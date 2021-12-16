#!/usr/bin/env python
# coding: utf-8

# In[3]:


# Loading necessary libraries
#from Bio import SeqIO
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
from sklearn import decomposition
from sklearn_som.som import SOM
from sklearn.manifold import Isomap
import geopandas as gpd
from sklearn.mixture import GaussianMixture
from scipy.special import logit
from sklearn.model_selection import train_test_split
from pykrige.rk import RegressionKriging
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import scale
from sklearn.metrics import mean_squared_error
from pykrige.ok import OrdinaryKriging
from sklearn.linear_model import LinearRegression


# In[60]:


# # Note - this blows up 4^k, so don't use this for long k-mers
# # Defining a function to generate the kmer frequency matrix for a genome with any given k (though it blows up exponentially so not too many)
# def calculate_kmer_frequencies(k,sequence):
#     kmer_library=itertools.product('ATCG',repeat=k)
#     output=[sequence.seq.count_overlap(''.join(kmer)) for kmer in kmer_library]
#     #for kmer in kmer_library:
#     #    kmer_count=sequence.seq.count_overlap(''.join(kmer))
#     #    output.append(kmer_count)
#     array_version=np.array(output)
#     return array_version


# # In[61]:


# ## DON'T RUN THIS I ALREADY SAVED THE OUTPUT
# # Here we are generating the tetramer matrix for all of the GOV viruses
# sequence_file=SeqIO.parse('GOV2_viral_populations_larger_than_5KB_or_circular.fasta','fasta')
# seq_names=[]
# output=[]
# tracker=0
# for record in sequence_file:
#     tracker+=1
#     seq_names.append(record.id)
#     kmers=calculate_kmer_frequencies(4,record)
#     output.append(kmers.transpose())
#     if tracker % 100 == 0:
#         print(str(tracker)+' sequences analyzed')
# output=np.stack(output,axis=0)
# print(output)


# # In[70]:


# ## ALSO DON'T RUN THIS UNLESS YOU RAN THE PREVIOUS CELL
# # Clearning up the output and writing out
# column_names=[''.join(item) for item in itertools.product('ATCG',repeat=4)]
# finished_frame=pd.DataFrame(output,index=seq_names,columns=column_names)
# finished_frame.to_csv('tetramer_frequencies.csv',index=True,header=True)
# finished_frame.index = finished_frame.index.astype('str')
# no_malaspina=finished_frame[finished_frame.index.str.contains('Station')]
# no_malaspina.to_csv('tetramer_gov_only.csv',index=True,header=True)


# In[4]:


# Reading in datasets that have already been generated
# genome_covariates.csv is part of code appendix A and includes the individual features for all genomes
# balanced_subsample.csv is part of code appendix D and is the data subset used to generate the regression kriging model
genome_data=pd.read_csv('/Users/rrodriguez77/Dropbox (GaTech)/viral_genome_stoichiometry/gov_2.0/genome_covariates.csv')
no_malaspina=pd.read_csv('/Users/rrodriguez77/Dropbox (GaTech)/viral_genome_stoichiometry/gov_2.0/tetramer_gov_only.csv',index_col=0)
normalized_frequencies=np.divide(no_malaspina.T,genome_data['length'].to_numpy()).T
#SOM_table = full_table.SOM
#full_table=example_genomes.set_index('Contig').join(SOM_table)
example_genomes=pd.read_csv('balanced_subsample.csv')
full_table=example_genomes.set_index('Contig').join(normalized_frequencies)
# Recalling the number of columns that go before the kmer frequencies
print(full_table.columns[-256:])
# Just making sure we don't miss any

# ROGELIO: TO AVOID PROBLEMS WHEN PLOTTING JUST CHANGE THE INDEX TO NUMERIC VALUES ALL INDEXES ARE NOW DIFFERENT WITH CONTIG WE HAVE DUPLICATE INDEXES
full_table.set_index(np.arange(20000), inplace = True)

# In[5]:


# Making a table that's numpy-array friendly to use for sklearn
kmers_only=full_table.iloc[:,-256:]
kmers_only
# Mkaing sure I know which columns are in there
full_table.columns[0:100]


# In[6]:


# Generating a GMM to delineate the bimodal distribution of C:N Ratios in the dataset
subpopulation_data=full_table[['c_to_n']].to_numpy()
# Learning the gaussian mixture
nitrogen_mixture=GaussianMixture(n_components=2,random_state=931).fit(subpopulation_data)
# Printing relevant metrics
print('BIC:')
print(nitrogen_mixture.bic(subpopulation_data))
print('GMM Means')
print(nitrogen_mixture.means_)
print('GMM Covariance matrix (variances)')
print(nitrogen_mixture.covariances_)
print('GMM Mixture Weights')
print(nitrogen_mixture.weights_)
# Assigning labels to the dataset
full_table['Distribution']=nitrogen_mixture.predict(full_table.loc[:,'c_to_n'].to_numpy().reshape(-1,1))
# Simulating a sample to compare to real data
points,labels=nitrogen_mixture.sample(20000)
# Plotting Results to Compare Generative Simulation to Observed Data
output_frame=pd.DataFrame(points,columns=['c_to_n'])
output_frame['Distribution']=labels
fig,ax=plt.subplots(1,2,sharey=True)
fig.set_size_inches(8,4)
ax[0].set_title('Simulated Gaussian Mixture')
ax[0].set_xlabel('C:N Ratio')
ax[1].set_title('Predicted Label for Real Data')
ax[1].set_xlabel('C:N Ratio')
sns.histplot(x=output_frame['c_to_n'],hue=output_frame['Distribution'],ax=ax[0],multiple='stack')
sns.histplot(x=full_table['c_to_n'],hue=full_table['Distribution'],ax=ax[1],multiple='stack')


# In[9]:


## Using the tetranucleotide frequencies to learn a map to show the differences in genome composition in a smaller
## dimensional space
kmer_som=SOM(m=2,n=2,dim=np.shape(kmers_only)[1],random_state=8131)
print('Learning Som')
kmer_som.fit(kmers_only.to_numpy())
print('Labeling with Som')
cluster_prediction=kmer_som.predict(kmers_only.to_numpy())
print('Learning Manifold with Isomap')
kmer_map=Isomap(n_components=2,p=1).fit_transform(kmers_only.to_numpy())
print('Plotting')
fig,ax=plt.subplots(2,1)
fig.set_size_inches(12,12)
ax[0].set_title('Isomap Reconstruction coded by SOM Cluster')
ax[1].set_title('C:N Ratios of SOM Clusters')
ax[1].set_xlabel('SOM Cluster ID')
sns.scatterplot(x=kmer_map[:,0],y=kmer_map[:,1],hue=cluster_prediction,palette='Set3',ax=ax[0])
sns.boxplot(x=cluster_prediction,y=full_table['c_to_n'],palette='Set3',ax=ax[1])
ax[1].set_ylabel('C:N Ratio')

sns.displot(x=cluster_prediction,hue=full_table['Distribution'],multiple='stack')
# Adding the SOM cluster labels to the dataframe
full_table['SOM_Cluster']=cluster_prediction

#full_table.to_csv('full_table.csv')


# In[155]:

full_table = pd.read_csv('full_table_4_SOMs.csv', index_col=0)

## Running regression again with the blocks of the cluster labels
## retrieving 2D locations
coords=pd.DataFrame()
coords['longitude']=full_table['longitude']+180
coords['latitude']=full_table['latitude']
sns.pairplot(coords)
#covar=full_table[['depth','nitrate', 'temperature', 'salinity']] # 
#covar=full_table[['depth', 'Distribution', 'SOM_Cluster']] # 
covar=full_table[['depth','nitrate', 'Distribution', 'SOM_Cluster', 'temperature', 'salinity']] # 
covar.insert(2,value=np.log(full_table['length']),column='log_length')
covar.insert(3,value=np.log(full_table['oxygen']),column='log_oxygen')
target=full_table['c_to_n']


# In[156]:


cov_train, cov_test, coord_train, coord_test, target_train, target_test = train_test_split(
    covar, coords, target, test_size=0.3, random_state=42
)
# This time because we have a mix of categorical variables and continuous variables with different bounds,
# instead of using a linear model we will use a random forest
regression_model=RandomForestRegressor(random_state=7645) 
lr_model = LinearRegression(normalize=True, copy_X = True, fit_intercept=True) 


# In[157]:


m_rk = RegressionKriging(regression_model = regression_model, n_closest_points = 10, coordinates_type = 'geographic', variogram_model = 'gaussian', nlags = 100)
m_rk.fit(cov_train.iloc[:,:].values, coord_train.iloc[:,:].values, target_train.iloc[:].values)
print("Regression Score (R^2): ", m_rk.regression_model.score(cov_test, target_test))
rk_score = m_rk.score(cov_test.iloc[:,:].values, coord_test.iloc[:,:].values, target_test.iloc[:].values)
print("RK score (R^2): ", rk_score)
print("RK Mean Squared Error: ", mean_squared_error(y_true = target_test.iloc[:].values, y_pred = m_rk.predict(cov_test.iloc[:,:].values, coord_test.iloc[:,:].values)))
predictions =  m_rk.predict(cov_test.iloc[:,:].values, coord_test.iloc[:,:].values)

# In[ ]:
# Final prediction of genomic C:N using SOM clusters as a new predictor variable

sns.set(rc = {'figure.figsize':(10,10)}, font_scale=1.5)
sns.set_style('whitegrid')
recovery=sns.scatterplot(y=predictions,x=target_test,hue=cov_test['SOM_Cluster'],palette='Set3')
recovery.set(xlabel='C:N Actual',ylabel='C:N Predicted')
plt.text(1.2, 1.375, r'$R^2$ = ' + "%.4f" % rk_score, fontsize = 14, weight = 'semibold')
plt.legend(loc='lower right', title = 'SOM')
#plt.show(recovery)

# In[158]:

# Plotting the variograms of regression residuals and the response variable

from pykrige.ok import OrdinaryKriging
from matplotlib.pyplot import figure

    
regression_prediction = m_rk.regression_model.predict(cov_train.iloc[:,:].values)
variogram_m = ['linear', 'power', 'gaussian', 'spherical', 'exponential']

data_norm = pd.read_csv('/Users/rrodriguez77/Dropbox (GaTech)/viral_genome_stoichiometry/gov_2.0/alldata_norm.csv', header=0)

figure(figsize=(12, 10))
for idx, m in enumerate(variogram_m):
    tmp = 231+idx
    
    # using full_table data
    # variogram of the regression residuals
    OK = OrdinaryKriging(coord_train.longitude, coord_train.latitude, target_train - regression_prediction, variogram_model = m, coordinates_type = 'geographic', nlags = 100)
    
    # variogram of the response variable
    #OK = OrdinaryKriging(coord_train.longitude, coord_train.latitude, target_train, variogram_model = m, coordinates_type = 'geographic', nlags = 100)
    
    #using Marian's normalized daya
    #OK = OrdinaryKriging(data_norm.longitude, data_norm.latitude, data_norm.total_n, variogram_model = m, coordinates_type = 'geographic', nlags = 70)
    variogram_data = OK.get_variogram_points();
    lags = variogram_data[0];
    semivariance_model = variogram_data[1];
    semivariance_data = OK.semivariance;
    plt.subplot(tmp)
    plt.plot(lags*111, semivariance_model, 'k-')
    plt.plot(lags*111, semivariance_data, 'ro')
    #plt.ylim(3.2e-5, 4.1e-5)
    plt.title('Variogram model: '+ m)
    plt.xlabel('distance (km)', fontsize= 12); plt.ylabel('semivariance',  fontsize= 12)




# In[259]:


m_rk.regression_model.get_params()


# In[262]:
unique_spots=coords.drop_duplicates()
spatial_residuals=m_rk.krige_residual(unique_spots.iloc[:,:].values)
xgrid=np.linspace(-180,180,360)+180
ygrid=np.linspace(-65,80,150)
full_grid=pd.DataFrame()
full_grid['lon']=np.repeat(xgrid,len(ygrid))
full_grid['lat']=np.tile(ygrid,len(xgrid))
grid_residuals=m_rk.krige_residual(full_grid.iloc[:,:].values)
spatial_residual_frame=pd.DataFrame()
grid_residuals_frame=pd.DataFrame()
spatial_residual_frame['longitude']=unique_spots['longitude']-180
spatial_residual_frame['latitude']=unique_spots['latitude']
spatial_residual_frame['value']=spatial_residuals
grid_residuals_frame['longitude']=full_grid['lon']-180
grid_residuals_frame['latitude']=full_grid['lat']
grid_residuals_frame['values']=grid_residuals

## Making maps ughhhhh
countries = gpd.read_file(
               gpd.datasets.get_path("naturalearth_lowres"))
countries = countries[(countries.name != "Antarctica") & (countries.name != "Fr. S. Antarctic Lands") & (countries.name!='Greenland')]
#countries=countries.to_crs(epsg=3832)
fig, ax = plt.subplots(figsize=(10,6))
ax.set_facecolor('white')
ax.set_title('Spatial Effect on Viral Genomic C:N')
grid_residuals_frame.plot(x='longitude',y='latitude',kind='scatter',ax=ax,c='values',cmap='viridis',s=0.45)
countries.plot(color="oldlace",edgecolor='oldlace',ax=ax,linewidth=0.4)
spatial_residual_frame.plot(x='longitude',y='latitude',kind='scatter',c='black',ax=ax)
#ax.set_xlim(left=-180,right=-20)
plt.show()



# In[ ]: Checking regression score on training and entire dataset
    
sample_data = pd.read_csv('full_table_4_SOMs.csv', index_col=0)

## retrieving 2D locations
coords=pd.DataFrame()
coords['longitude']=sample_data['longitude']+180
coords['latitude']=sample_data['latitude']
sns.pairplot(coords)
#covar=sample_data[['depth','nitrate']]
covar = sample_data[['depth','nitrate', 'salinity', 'temperature']] # 
covar.insert(2, value=np.log(sample_data['length']),column='log_length')
covar.insert(3, value=np.log(sample_data['oxygen']),column='log_oxygen')
#target=logit(sample_data['GC']/100)
target = sample_data['c_to_n']
#covar = scale(covar)



cov_train, cov_test, coord_train, coord_test, target_train, target_test = train_test_split(
    covar, coords, target, test_size=0.3, random_state=42
)

# checking how good is regression fit on training and entire dataset
lr_model = LinearRegression(normalize = True, copy_X = True, fit_intercept = True)
lr_model.fit(cov_train, target_train)
print("Regression Score of trainig set: ", lr_model.score(cov_train, target_train))
lr_model.fit(covar, target)
print("Regression Score entire dataset: ", lr_model.score(covar, target))

# checking how good is regression fit on training and entire dataset using SOM clusters
covar.insert(4, value= sample_data['SOM_Cluster'],column='SOM_Cluster')
covar.insert(5, value= sample_data['Distribution'],column='Distribution')
cov_train, cov_test, coord_train, coord_test, target_train, target_test = train_test_split(
    covar, coords, target, test_size=0.3, random_state=42
)
#lr_model = LinearRegression(normalize = True, copy_X = True, fit_intercept = True)
lr_model = RandomForestRegressor(random_state=7645)
lr_model.fit(cov_train, target_train)
print("Regression Score of trainig set: ", lr_model.score(cov_train, target_train))
lr_model.fit(covar, target)
print("Regression Score entire dataset: ", lr_model.score(covar, target))




# In[14]:


#Now that we know both of these designations (SOM cluster and C:N GMM label) are ecologically significant
#we can start investigating the details of what biogeographies characterize each group
# First step, let's assign labels to ALL the genomes >:D
all_som_assignments=kmer_som.predict(normalized_frequencies.to_numpy())
all_gmm_assignments_prob=nitrogen_mixture.predict_proba(genome_data['c_to_n'].to_numpy().reshape(-1,1))
all_gmm_assignments_binary=nitrogen_mixture.predict(genome_data['c_to_n'].to_numpy().reshape(-1,1))


# In[46]:


#cluster_labels_frame=pd.DataFrame()
#cluster_labels_frame['Contig']=genome_data['ID']
#cluster_labels_frame['SOM']=all_som_assignments
#cluster_labels_frame['GMM_decision']=all_gmm_assignments_binary
#cluster_labels_frame['GMM_highN']=all_gmm_assignments_prob[:,0]
#cluster_labels_frame['GMM_lowN']=all_gmm_assignments_prob[:,1]
sns.histplot(data=cluster_labels_frame,x='SOM',hue='GMM_decision',multiple='stack')
sns.jointplot(data=cluster_labels_frame,x='GMM_highN',y='GMM_lowN',kind='hex')
# Relatively, there are extremely few ambiguous cases in the dataset


# In[99]:


# Now relating that to the abundance data in the different samples
#abundance_data=pd.read_csv('alldata.csv')
#abundance_data_subset=abundance_data[['Contig','length','GC','c_to_n','total_n','abundance','level_code','station_id','latitude','longitude','depth']]
#abundance_and_cluster=abundance_data_subset.set_index('Contig').join(cluster_labels_frame.set_index('Contig'))
abundance_and_cluster['level_code'].value_counts()
#abundance_and_cluster['SOM'] = pd.Categorical(abundance_and_cluster['SOM'], np.arange(12)).as_ordered()


# In[114]:


# Looking at differences with respect to depth layers
sns.histplot(data=abundance_and_cluster,y='SOM',weights='abundance',hue='GMM_decision',multiple='stack')


# In[148]:


## Plotting the distribution of surface viruses
srf=sns.displot(data=abundance_and_cluster[abundance_and_cluster['level_code'].isin(['SUR'])],
            height=6,
            hue='SOM',
            y='latitude',
            weights='abundance',
            bins=50,
            palette='Set3',
            multiple='stack')
srf.fig.suptitle('Surface Viruses')
srf.ax.set_xlabel('Abundance (Read Counts)')
srf.ax.set_ylabel('Latitude')


# In[179]:


## Plotting the distribution of surface viruses
dcm=sns.displot(data=abundance_and_cluster[abundance_and_cluster['level_code'].isin(['DCM'])],
            height=6,
            hue='SOM',
            y='latitude',
            weights='abundance',
            bins=50,
            palette='Set3',
            multiple='stack')
dcm.fig.suptitle('DCM Viruses')
dcm.ax.set_xlabel('Abundance (Read Counts)')
dcm.ax.set_ylabel('Latitude')


# In[186]:


## Plotting the distribution of surface viruses
mes=sns.displot(data=abundance_and_cluster[abundance_and_cluster['level_code'].isin(['MES'])],
            height=6,
            hue='SOM',
            y='latitude',
            weights='abundance',
            bins=50,
            palette='Set3',
            multiple='stack')
mes.fig.suptitle('Mesopelagic Viruses')
mes.ax.set_xlabel('Abundance (Read Counts)')
mes.ax.set_ylabel('Latitude')


# In[170]:


mean_GMM=abundance_and_cluster.groupby(['station_id']).apply(lambda x: (x['GMM_decision']*x['abundance']).sum()/x['abundance'].sum())
unique_cluster_spots=abundance_and_cluster.drop_duplicates('station_id').set_index('station_id')
unique_cluster_spots['mean_GMM']=mean_GMM


# In[178]:


## Making maps ughhhhh
countries = gpd.read_file(
               gpd.datasets.get_path("naturalearth_lowres"))
countries = countries[(countries.name != "Antarctica") & (countries.name != "Fr. S. Antarctic Lands") & (countries.name!='Greenland')]
#countries=countries.to_crs(epsg=3832)
fig, ax = plt.subplots(figsize=(10,6))
ax.set_facecolor('white')
ax.set_title('Spatial Effect on Viral Genomic C:N')
grid_residuals_frame.plot(x='longitude',y='latitude',kind='scatter',ax=ax,c='values',cmap='viridis',s=0.45)
countries.plot(color="oldlace",edgecolor='oldlace',ax=ax,linewidth=0.4)
unique_cluster_spots.plot(x='longitude',y='latitude',kind='scatter',c='mean_GMM',ax=ax)
plt.show()


# In[ ]:




