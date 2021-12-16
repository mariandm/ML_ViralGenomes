#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 19:42:22 2021

@author: rrodriguez77
"""

from sklearn.cluster import KMeans
import numpy as np
from matplotlib.pyplot import figure
import pandas as pd
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


data_norm = pd.read_csv('/Users/rrodriguez77/Dropbox (GaTech)/viral_genome_stoichiometry/gov_2.0/alldata_norm.csv', header=0)


kmeans = KMeans(n_clusters=4, random_state=0).fit(data_norm[['c_to_n', 'GC']])
kmeans.labels_


# Creating dataset
y = data_norm.depth
x = data_norm.longitude
z = data_norm.latitude
#size= np.array(data_norm.GC - np.mean(data_norm.GC))/np.std(data_norm.GC)
size= np.array((1/data_norm.c_to_n) - np.mean(1/data_norm.c_to_n)) / np.std(1/data_norm.c_to_n)
size = size + abs(min(size)) + 0.1

# Creating figure
fig = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d")
# Creating plot
sc = ax.scatter3D(x, y, z, c = kmeans.labels_, cmap = 'Set2', s = 20*size, alpha=1)
plt.title("Stations clustering")
plt.xlabel('longitude'); plt.ylabel('depth');
ax.set_zlabel('latitude')
plt.legend(*sc.legend_elements(), bbox_to_anchor=(0.9, 1.1), loc=2)
# show plot
plt.show()

depth_level = [data_norm.station_abundance[i].split("_")[1] for i in range(data_norm.shape[0])]
indices_sur = [i for i, x in enumerate(depth_level) if x == "SUR"]
indices_dcm = [i for i, x in enumerate(depth_level) if x == "DCM"]
indices_mes = [i for i, x in enumerate(depth_level) if x == "MES"]
indices_mix = [i for i, x in enumerate(depth_level) if x == "MIX"]
indices_zzz = [i for i, x in enumerate(depth_level) if x == "ZZZ"]



import matplotlib.pyplot as plt
import seaborn as sns
import geopandas as gpd
countries = gpd.read_file(
              gpd.datasets.get_path("naturalearth_lowres"))
countries = countries[(countries.name != "Antarctica") & (countries.name != "Fr. S. Antarctic Lands") & (countries.name!='Greenland')]


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
ax1.set_facecolor('white')
ax1.set_title('SUR')
countries.plot(color="oldlace",edgeincolor='oldlace',ax=ax1,linewidth=1)
ax1.scatter(x[indices_sur], z[indices_sur], c = kmeans.labels_[indices_sur], cmap = 'Set2', s = 20*size[indices_sur])
ax1.set_xlabel('Longitude', fontsize = 13)
ax1.set_ylabel('Latitude', fontsize = 13)
ax1.text(-180, -55, 'Depth: ' + '%.1f' % min(y[indices_sur]) + ' to ' + '%.1f' % max(y[indices_sur]) + ' m' , fontsize = 12, fontweight = 'bold')

ax2.set_facecolor('white')
ax2.set_title('DCM')
countries.plot(color="oldlace",edgecolor='oldlace',ax=ax2,linewidth=1)
ax2.scatter(x[indices_dcm], z[indices_dcm], c = kmeans.labels_[indices_dcm], cmap = 'Set2', s = 20*size[indices_dcm])
ax2.set_xlabel('Longitude', fontsize = 13)
ax2.set_ylabel('Latitude', fontsize = 13)
ax2.text(-180, -55, 'Depth: ' + '%.1f' % min(y[indices_dcm]) + ' to ' + '%.1f' % max(y[indices_dcm]) + ' m' , fontsize = 12, fontweight = 'bold')


ax3.set_facecolor('white')
ax3.set_title('MIX')
countries.plot(color="oldlace",edgecolor='oldlace',ax=ax3,linewidth=1)
ax3.scatter(x[indices_mix], z[indices_mix], c = kmeans.labels_[indices_mix], cmap = 'Set2', s = 20*size[indices_mix])
ax3.set_xlabel('Longitude', fontsize = 13)
ax3.set_ylabel('Latitude', fontsize = 13)
ax3.text(-180, -55, 'Depth: ' + '%.1f' % min(y[indices_mix]) + ' to ' + '%.1f' % max(y[indices_mix]) + ' m' , fontsize = 12, fontweight = 'bold')


ax4.set_facecolor('white')
ax4.set_title('MES')
countries.plot(color="oldlace",edgecolor='oldlace',ax=ax4,linewidth=1)
sc = ax4.scatter(x[indices_mes], z[indices_mes], c = kmeans.labels_[indices_mes], cmap = 'Set2', s = 20*size[indices_mes])
ax4.set_xlabel('Longitude', fontsize = 13)
ax4.set_ylabel('Latitude', fontsize = 13)
ax4.text(-180, -55, 'Depth: ' + '%.1f' % min(y[indices_mes]) + ' to ' + '%.1f' % max(y[indices_mes]) + ' m' , fontsize = 12, fontweight = 'bold')
plt.legend(*sc.legend_elements(), bbox_to_anchor=(1, 1.1), loc=2, title = 'Clusters', fontsize = 13)
plt.show()


## SOM clustering marian
station_clust_som = pd.read_csv('stationgroup.csv', header=0)
station_clust_som.dropna(subset = ["latitude"], inplace=True)



fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
ax1.set_facecolor('white')
ax1.set_title('SUR')
countries.plot(color="oldlace",edgecolor='oldlace',ax=ax1,linewidth=1)
ax1.scatter(x[indices_sur], z[indices_sur], c = station_clust_som.kmeans.iloc[indices_sur], cmap = 'Accent', s = 20*size[indices_sur])
ax1.set_xlabel('Longitude', fontsize = 13)
ax1.set_ylabel('Latitude', fontsize = 13)
ax1.text(-180, -55, 'Depth: ' + '%.1f' % min(y[indices_sur]) + ' to ' + '%.1f' % max(y[indices_sur]) + ' m' , fontsize = 12, fontweight = 'bold')

ax2.set_facecolor('white')
ax2.set_title('DCM')
countries.plot(color="oldlace",edgecolor='oldlace',ax=ax2,linewidth=1)
ax2.scatter(x[indices_dcm], z[indices_dcm], c =station_clust_som.kmeans.iloc[indices_dcm], cmap = 'Accent', s = 20*size[indices_dcm])
ax2.set_xlabel('Longitude', fontsize = 13)
ax2.set_ylabel('Latitude', fontsize = 13)
ax2.text(-180, -55, 'Depth: ' + '%.1f' % min(y[indices_dcm]) + ' to ' + '%.1f' % max(y[indices_dcm]) + ' m' , fontsize = 12, fontweight = 'bold')


ax3.set_facecolor('white')
ax3.set_title('MIX')
countries.plot(color="oldlace",edgecolor='oldlace',ax=ax3,linewidth=1)
ax3.scatter(x[indices_mix], z[indices_mix],  c = ['#f0027f', '#fdc086', '#7fc97f', '#fdc086', '#f0027f' ], s = 20*size[indices_mix])## ['#7fc97f', '#fdc086', '#f0027f', '#666666']
ax3.set_xlabel('Longitude', fontsize = 13)
ax3.set_ylabel('Latitude', fontsize = 13)
ax3.text(-180, -55, 'Depth: ' + '%.1f' % min(y[indices_mix]) + ' to ' + '%.1f' % max(y[indices_mix]) + ' m' , fontsize = 12, fontweight = 'bold')


ax4.set_facecolor('white')
ax4.set_title('MES')
countries.plot(color="oldlace",edgecolor='oldlace',ax=ax4,linewidth=1)
sc = ax4.scatter(x[indices_mes], z[indices_mes], c = station_clust_som.kmeans.iloc[indices_mes], cmap = 'Accent', s = 20*size[indices_mes])
ax4.set_xlabel('Longitude', fontsize = 13)
ax4.set_ylabel('Latitude', fontsize = 13)
ax4.text(-180, -55, 'Depth: ' + '%.1f' % min(y[indices_mes]) + ' to ' + '%.1f' % max(y[indices_mes]) + ' m' , fontsize = 12, fontweight = 'bold')
plt.legend(*sc.legend_elements(), bbox_to_anchor=(1, 1.1), loc=2, title = 'Clusters', fontsize = 13)
plt.show()


# In[ ]: Stations clustering map with environmental nitrate concentration levels

from pykrige.ok import OrdinaryKriging
from matplotlib.pyplot import figure

    
# variogram_m = ['linear', 'power', 'gaussian', 'spherical', 'exponential']


# figure(figsize=(12, 10))
# for idx, m in enumerate(variogram_m):
#     tmp = 231+idx
    
#     # using full_table data
#     #OK = OrdinaryKriging(coord_train.longitude, coord_train.latitude, target_train - regression_prediction, variogram_model = m, coordinates_type = 'geographic', nlags = 100)
#     OK = OrdinaryKriging(coords.longitude, coords.latitude, covar.nitrate, variogram_model = m, coordinates_type = 'geographic', nlags = 70)
    
#     variogram_data = OK.get_variogram_points();
#     lags = variogram_data[0];
#     semivariance_model = variogram_data[1];
#     semivariance_data = OK.semivariance;
#     plt.subplot(tmp)
#     plt.plot(lags*111, semivariance_model, 'k-')
#     plt.plot(lags*111, semivariance_data, 'ro')
#     #plt.ylim(3.2e-5, 4.1e-5)
#     plt.title('Variogram model: '+ m)
#     plt.xlabel('distance (km)', fontsize= 12); plt.ylabel('semivariance',  fontsize= 12)


sample_data = pd.read_csv('full_table_4_SOMs.csv', index_col=0)

    
# Separate data by depth level
index_sur = np.where(sample_data.level_code == 'SUR')
index_dcm = np.where(sample_data.level_code == 'DCM')
index_mix = np.where(sample_data.level_code == 'MIX')
index_mes = np.where(sample_data.level_code == 'MES')

data_sur = sample_data.iloc[index_sur[0],:]   
data_dcm = sample_data.iloc[index_dcm[0],:] 
data_mix = sample_data.iloc[index_mix[0],:] 
data_mes = sample_data.iloc[index_mes[0],:]


xgrid=np.linspace(-180,180, 100)+180
ygrid=np.linspace(-65,80,100)

OK = OrdinaryKriging(data_sur.longitude, data_sur.latitude, data_sur.nitrate, variogram_model = 'spherical', coordinates_type = 'geographic', nlags = 18, enable_plotting = True)
nitrate_val_sur, ss = OK.execute('grid', xgrid, ygrid)

OK = OrdinaryKriging(data_dcm.longitude, data_dcm.latitude, data_dcm.nitrate, variogram_model = 'spherical', coordinates_type = 'geographic', nlags = 18, enable_plotting = True)
nitrate_val_dcm, ss = OK.execute('grid', xgrid, ygrid)

OK = OrdinaryKriging(data_mix.longitude, data_mix.latitude, data_mix.nitrate, variogram_model = 'spherical', coordinates_type = 'geographic', nlags = 18, enable_plotting = True)
nitrate_val_mix, ss = OK.execute('grid', xgrid, ygrid)

OK = OrdinaryKriging(data_mes.longitude, data_mes.latitude, data_mes.nitrate, variogram_model = 'spherical', coordinates_type = 'geographic', nlags = 18, enable_plotting = True)
nitrate_val_mes, ss = OK.execute('grid', xgrid, ygrid)

mes_nitrate = np.array(data_mes.nitrate)
sur_nitrate = np.array(data_sur.nitrate)
mix_nitrate = np.array(data_mix.nitrate)
dcm_nitrate = np.array(data_dcm.nitrate)

data_mes[['longitude','latitude']].iloc[np.where(mes_nitrate > 17)].drop_duplicates()

import matplotlib.pyplot as plt
import seaborn as sns
import geopandas as gpd
countries = gpd.read_file(
              gpd.datasets.get_path("naturalearth_lowres"))
countries = countries[(countries.name != "Antarctica") & (countries.name != "Fr. S. Antarctic Lands") & (countries.name!='Greenland')]


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
ax1.set_facecolor('white')
ax1.set_title('SUR')
sur = ax1.imshow(nitrate_val_sur, extent= [-180, 180, -65, 80], origin = 'lower')
fig.colorbar(sur, ax=ax1)
ax1.set_ylim(-65, 80)
countries.plot(color="oldlace",edgecolor='oldlace',ax=ax1,linewidth=1)
ax1.scatter(x[indices_sur], z[indices_sur], c = station_clust_som.kmeans.iloc[indices_sur], cmap = 'Accent', s = 20*size[indices_sur])
ax1.set_xlabel('Longitude', fontsize = 13)
ax1.set_ylabel('Latitude', fontsize = 13)
ax1.text(-180, -55, 'Depth: ' + '%.1f' % min(y[indices_sur]) + ' to ' + '%.1f' % max(y[indices_sur]) + ' m' , fontsize = 12, fontweight = 'bold')

ax2.set_facecolor('white')
ax2.set_title('DCM')
dcm = ax2.imshow(nitrate_val_dcm, extent= [-180, 180, -65, 80], origin = 'lower')
fig.colorbar(dcm, ax=ax2)
ax2.set_ylim(-65, 80)
countries.plot(color="oldlace",edgecolor='oldlace',ax=ax2,linewidth=1)
ax2.scatter(x[indices_dcm], z[indices_dcm], c =station_clust_som.kmeans.iloc[indices_dcm], cmap = 'Accent', s = 20*size[indices_dcm])
ax2.set_xlabel('Longitude', fontsize = 13)
ax2.set_ylabel('Latitude', fontsize = 13)
ax2.text(-180, -55, 'Depth: ' + '%.1f' % min(y[indices_dcm]) + ' to ' + '%.1f' % max(y[indices_dcm]) + ' m' , fontsize = 12, fontweight = 'bold')


ax3.set_facecolor('white')
ax3.set_title('MIX')
mix = ax3.imshow(nitrate_val_mix, extent= [-180, 180, -65, 80], origin = 'lower')
fig.colorbar(mix, ax=ax3)
ax3.set_ylim(-65, 80)
countries.plot(color="oldlace",edgecolor='oldlace',ax=ax3,linewidth=1)
ax3.scatter(x[indices_mix], z[indices_mix],  c = ['#f0027f', '#fdc086', '#7fc97f', '#fdc086', '#f0027f' ], s = 20*size[indices_mix])## ['#7fc97f', '#fdc086', '#f0027f', '#666666']
ax3.set_xlabel('Longitude', fontsize = 13)
ax3.set_ylabel('Latitude', fontsize = 13)
ax3.text(-180, -55, 'Depth: ' + '%.1f' % min(y[indices_mix]) + ' to ' + '%.1f' % max(y[indices_mix]) + ' m' , fontsize = 12, fontweight = 'bold')


ax4.set_facecolor('white')
ax4.set_title('MES')
mes = ax4.imshow(nitrate_val_mes, extent= [-180, 180, -65, 80], origin = 'lower')
fig.colorbar(mes, ax=ax4)
ax4.set_ylim(-65, 80)
countries.plot(color="oldlace",edgecolor='oldlace',ax=ax4,linewidth=1)
sc = ax4.scatter(x[indices_mes], z[indices_mes], c = station_clust_som.kmeans.iloc[indices_mes], cmap = 'Accent', s = 20*size[indices_mes])
ax4.set_xlabel('Longitude', fontsize = 13)
ax4.set_ylabel('Latitude', fontsize = 13)
ax4.text(-180, -55, 'Depth: ' + '%.1f' % min(y[indices_mes]) + ' to ' + '%.1f' % max(y[indices_mes]) + ' m' , fontsize = 12, fontweight = 'bold')
plt.legend(*sc.legend_elements(), bbox_to_anchor=(1.2, 0.8), loc=2, title = 'Clusters', fontsize = 13)
plt.show()








