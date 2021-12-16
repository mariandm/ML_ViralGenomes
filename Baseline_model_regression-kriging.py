#!/usr/bin/env python
# coding: utf-8

# In[1]:


## Setting up a minimal example for doing krigging regression
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import random
from sklearn.linear_model import LinearRegression
from sklearn.neural_network import MLPRegressor
from sklearn.linear_model import LogisticRegression
from scipy.special import logit
from scipy.special import expit
from sklearn.model_selection import train_test_split
from pykrige.rk import RegressionKriging
from sklearn.ensemble import RandomForestRegressor



# In[2]:


## Reading in viral genome abundance and distribution data
data_file = pd.read_csv('/Users/rrodriguez77/Dropbox (GaTech)/viral_genome_stoichiometry/gov_2.0/alldata_111821.csv', header=0)


# In[3]:


all_stations=data_file[['latitude','depth']].drop_duplicates(ignore_index=True)
all_stations=all_stations[all_stations['latitude'].isna()==False]
#station_sequences=data_file[(data_file['latitude']==all_stations['latitude'][0]) & (data_file['depth']==all_stations['depth'][0])]
#station_sequences=station_sequences.assign(sample_probability=station_sequences['abundance']/sum(station_sequences['abundance']))
#n_seq=min(1000,np.shape(station_sequences)[0])
#sns.histplot(station_sequences['sample_probability'])
#sampled_subset=np.random.choice(np.arange(np.shape(station_sequences)[0]),n_seq, replace=False, p=station_sequences['sample_probability'])
#output=station_sequences.loc[sampled_subset,:]
output = pd.DataFrame()
for index,row in all_stations.iterrows():
    station_data=data_file[(data_file['latitude']==row['latitude']) & (data_file['depth']==row['depth'])]
    station_data=station_data.assign(sample_probability=station_data['abundance']/sum(station_data['abundance']))
    n_seq=min(1000,np.shape(station_data)[0])
    subsamp= np.random.choice(np.arange(np.shape(station_data)[0]), n_seq, replace=False, p=station_data['sample_probability'])
    output= output.append(station_data.iloc[subsamp,:],ignore_index=True)
print(output)
sns.pairplot(output[['latitude','longitude','depth']])


# In[4]:


## Taking a small subset to just test the syntax
random.seed(980531)
## Get a line in to make sure no negative nitrate values
# Using all data smaller_data=data_file[data_file['nitrate'].notna()]
# Using blocked subsampled data
smaller_data=output[output['nitrate'].notna()]
smaller_data=smaller_data.sample(frac=1).reset_index(drop=True)
sns.pairplot(smaller_data[['latitude','longitude','depth']])
print(np.shape(smaller_data))
sample_rows=random.sample(range(len(smaller_data)),20000)
sample_data=smaller_data.iloc[sample_rows,:]
sample_data.head()
sns.pairplot(sample_data[['latitude','longitude','depth']])
#sns.pairplot(sample_data[['length','GC','c_to_n','abundance','latitude','longitude','depth','darwin_nitrate']])


# In[5]:

# Selection of predictor variables and target variable
    
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

# In[6]:

# Split into training and testing sets
cov_train, cov_test, coord_train, coord_test, target_train, target_test = train_test_split(
    covar, coords, target, test_size=0.3, random_state=42
)
lr_model = LinearRegression(normalize = True, copy_X = True, fit_intercept = True)
#rf_model=RandomForestRegressor(n_estimators = 200) 

target_train.iloc[:].values


# In[7]:

# Baseline model evaluation using linear regression
m_rk = RegressionKriging(regression_model=lr_model, n_closest_points=10, coordinates_type='geographic',variogram_model='gaussian', nlags = 70)
m_rk.fit(cov_train.iloc[:,:].values, coord_train.iloc[:,:].values, target_train.iloc[:].values)
print("Regression Score: ", m_rk.regression_model.score(cov_test, target_test))
print("RK score: ", m_rk.score(cov_test.iloc[:,:].values, coord_test.iloc[:,:].values, target_test.iloc[:].values))


# In[8]:


m_rk.regression_model.coef_


# In[9]:


m_rk.regression_model.intercept_


# In[10]:


unique_spots=coords.drop_duplicates()
spatial_residuals=m_rk.krige_residual(unique_spots.iloc[:,:].values)
xgrid=np.linspace(-180,180,360)+180
ygrid=np.linspace(-65,80,150)
full_grid=pd.DataFrame()
full_grid['lon']=np.repeat(xgrid,len(ygrid))
full_grid['lat']=np.tile(ygrid,len(xgrid))
grid_residuals=m_rk.krige_residual(full_grid.iloc[:,:].values)


# In[11]:


spatial_residual_frame=pd.DataFrame()
grid_residuals_frame=pd.DataFrame()
spatial_residual_frame['longitude']=unique_spots['longitude']-180
spatial_residual_frame['latitude']=unique_spots['latitude']
spatial_residual_frame['value']=(spatial_residuals)
grid_residuals_frame['longitude']=full_grid['lon']-180
grid_residuals_frame['latitude']=full_grid['lat']
grid_residuals_frame['values']=(grid_residuals)


# In[12]:

# Spatial effect on genomic C:N prediction
import matplotlib.pyplot as plt
import seaborn as sns
import geopandas as gpd
countries = gpd.read_file(
              gpd.datasets.get_path("naturalearth_lowres"))
countries = countries[(countries.name != "Antarctica") & (countries.name != "Fr. S. Antarctic Lands") & (countries.name!='Greenland')]
#countries=countries.to_crs(epsg=3832)
fig, ax = plt.subplots(figsize=(10,6))
ax.set_facecolor('white')
ax.set_title('Spatial Effect on Viral C:N ratio')
grid_residuals_frame.plot(x='longitude',y='latitude',kind='scatter',ax=ax, c='values',cmap='viridis',s=0.45, vmin=0,vmax=0.045) # vmin=0.4,vmax=0.6
countries.plot(color="oldlace",edgecolor='oldlace',ax= ax,linewidth=0.4)
spatial_residual_frame.plot(x='longitude',y='latitude',kind='scatter',ax= ax, c='black')
#ax.set_xlim(left=-180,right=-20)
plt.show()



# In[ ]:

# In this section we evaluate different regression methods used
# by the Regression-kriging model to predict GC content.
# These are the regression methods we evaluate: Support Vector Regression, 
# Random Forest, K-nearest neighbors regression, and Linear regression
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import scale
from sklearn.metrics import mean_squared_error

# Support vector regression
svr_model = SVR(C = 1, gamma = "scale")
# Random forest
rf_model = RandomForestRegressor(n_estimators = 200) 
# K-nearest neighbors
kneigh_model = KNeighborsRegressor(n_neighbors = 10, weights = 'distance') 
# Linear regression model
lr_model = LinearRegression(normalize=True, copy_X = True, fit_intercept=True) 

models = [svr_model, rf_model, kneigh_model, lr_model]


prediction_models = pd.DataFrame()
RK_r2 = {}

for m in models:
    print("=" * 40)
    print("regression model:", m.__class__.__name__)
    m_rk = RegressionKriging(regression_model = m, n_closest_points = 10, coordinates_type = 'geographic', variogram_model = 'gaussian', nlags = 70)
    m_rk.fit(cov_train.iloc[:,:].values, coord_train.iloc[:,:].values, target_train.iloc[:].values)
    print("Regression Score (R^2): ", m_rk.regression_model.score(cov_test, target_test))
    rk_score = m_rk.score(cov_test.iloc[:,:].values, coord_test.iloc[:,:].values, target_test.iloc[:].values)
    print("RK score (R^2): ", rk_score)
    print("RK Mean Squared Error: ", mean_squared_error(y_true = target_test.iloc[:].values, y_pred = m_rk.predict(cov_test.iloc[:,:].values, coord_test.iloc[:,:].values)))
    prediction_models[m.__class__.__name__] =  m_rk.predict(cov_test.iloc[:,:].values, coord_test.iloc[:,:].values)
    RK_r2[m.__class__.__name__] = rk_score


from matplotlib import pyplot as plt
from matplotlib.pyplot import figure

# Plot True vs Predicted GC content. For the different regression
# methods used above. Also, show the coefficient of determination
#  of the prediction
figure(figsize=(12, 12))
plt.subplot(221)
plt.scatter(target_test.iloc[:].values, prediction_models['LinearRegression'], c = 'k')
plt.plot(np.linspace(0,2,100), np.linspace(0,2,100), 'k--')
plt.title('Linear regression')
plt.xlim(1.05,1.49); plt.ylim(1.05,1.49)
plt.xlabel('True C:N', fontsize= 12); plt.ylabel('Predicted C:N',  fontsize= 12)
plt.text(1.07, 1.4, r'$R^2$ = ' + "%.4f" % RK_r2['LinearRegression'], fontsize = 12)


plt.subplot(222)
plt.scatter(target_test.iloc[:].values, prediction_models['SVR'], c = 'k')
plt.plot(np.linspace(0,2,100), np.linspace(0,2,100), 'k--')
plt.title('SVR')
plt.xlim(1.05,1.49); plt.ylim(1.05,1.49)
plt.xlabel('True C:N', fontsize= 12); plt.ylabel('Predicted C:N',  fontsize= 12)
plt.text(1.07, 1.4, r'$R^2$ = ' + "%.4f" % RK_r2['SVR'], fontsize = 12)


plt.subplot(223)
plt.scatter(target_test.iloc[:].values, prediction_models['RandomForestRegressor'], c = 'k')
plt.plot(np.linspace(0,2,100), np.linspace(0,2,100), 'k--')
plt.title('Random Forest')
plt.xlim(1.05,1.49); plt.ylim(1.05,1.49)
plt.xlabel('True C:N', fontsize= 12); plt.ylabel('Predicted C:N',  fontsize= 12)
plt.text(1.07, 1.4, r'$R^2$ = ' + "%.4f" % RK_r2['RandomForestRegressor'], fontsize = 12)


plt.subplot(224)
plt.scatter(target_test.iloc[:].values, prediction_models['KNeighborsRegressor'], c = 'k')
plt.plot(np.linspace(0,2,100), np.linspace(0,2,100), 'k--')
plt.title('K-nearest neighbors')
plt.xlim(1.05,1.49); plt.ylim(1.05,1.49)
plt.xlabel('True C:N', fontsize= 12); plt.ylabel('Predicted C:N',  fontsize= 12)
plt.text(1.07, 1.4, r'$R^2$ = ' + "%.4f" % RK_r2['KNeighborsRegressor'], fontsize = 12)


#plt.savefig('linkage_comparison.eps')

