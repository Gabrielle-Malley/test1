# Read Me

#Module Name: Grassland Curing Data Time Series

#Date of creation: 17/05/2022

#Name of creator of module : GM

#History of modification:
    #Last modification: 19/05
    #Modifications: Raw script to Juypter notebook 
    #- added a loop for multiple fwd and changes a few variable names
    #- added the curing data from years 2000-2019 and plotted
    #  

#Summary of what the module does: 

#Reads fire weather district (fwd) shape files as a dataframe 

#Reads in avaliable MODIS curing data from 2015-2020

#Loops over the day of year, for each year every 8 days

#MODIS data opened using xarray (xr), mean curing values are determined for each fire weather district

#Generates a basic plot for each fire weather district from 2002-2019 (not beginning at 2000 as the years 2000-2001 are incomplete

#typical calls/keypaths: 

#curing_file_path = '/g/data/er8/global_challenge/curing'

#fwd_shapefile_path = '/g/data/ct18/dw4060/viirs/inputfiles/boundaries/FWFDistricts_MapID_WGS84_May2022.shp'

#%%
#import libraries 
import os, sys
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import geopandas as gpd
import xarray as xr
import datetime
import glob
import os


#read in data
#read in fire weather district (fwd) shape files 
fwd_shapefile_path = '/g/data/ct18/dw4060/viirs/inputfiles/boundaries/FWFDistricts_MapID_WGS84_May2022.shp'
df = gpd.read_file(fwd_shapefile_path)

curing_file_path = '/g/data/er8/global_challenge/curing'

#define fire weather districts and make empty lists
fwd_districts = ["Wimmera", "Mallee", "South West", "Northern Country", "North Central", "Central", "North East", "West and South Gippsland", "East Gippsland"]
total_data = []
time_array = []
mean_curing_value = []

#get data for plots
for name in range(len(fwd_districts)):
    fwd_geometry = df[df['DIST_NAME'] == fwd_districts[name]].geometry.values
    year = 2002
    t0 = datetime.datetime(year, 4, 6)
    time_array = []
    mean_curing_value = []

    for doy in range(0, 6955, 7):
        current_date = t0 + datetime.timedelta(days=doy)
        #curing_modis-mapvictoria_8day_500m_vic_%s.tif
        fns = glob.glob(os.path.join(curing_file_path, "curing_modis-mapvictoria_8day_500m_vic_%s.tif" % (current_date.strftime('%Y%m%d')+'-'+((current_date + datetime.timedelta(7)).strftime('%Y%m%d')))))
        print(curing_file_path)
        print(current_date.strftime('%Y%m%d')+'-'+((current_date + datetime.timedelta(6)).strftime('%Y%m%d')))
        #curing-filled_modis-mapvictoria_8day_500m_vic_%s.ti
        #print("fns", fns)
    #cheking file exists
        if len(fns) == 1: 
        #opening the file using xr
            ds = xr.open_dataset(fns[0])
            #taking the grassland curing data (gci) and clipping I with the geometry of the fwds
            gci_clipped_to_fwd = ds.rio.clip(fwd_geometry, df.crs)
            #calculte the mean for the clipped gci
            #determine mean curing values for the all the fwds and structure them as an array starting from the first idex
            mean_curing_value.append(gci_clipped_to_fwd.mean().to_array().values[0])
            #add current date variable to the time arrazy variable 
            time_array.append(current_date)
            
#plotting out array (you can also stor you plot as a variable/object to manipulate it further)
    plt.plot(time_array,mean_curing_value)
    plt.title(fwd_districts[name])
    plt.xlabel('Time (Year-Month)')
    plt.ylabel('Mean Grassland Curing Value (%)')
    plt.show()
    total_data.append([time_array, mean_curing_value])
#plt.savefig("curing_timeline.png")

print(total_data)
    
