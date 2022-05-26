
import os, sys
import xarray as xr
import numpy as np
#import numpy.ma as ma
import pandas as pd
import glob
#from zipfile import ZipFile
from osgeo import gdal
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import rioxarray as rio 
import rasterio
from rasterio.plot import show
import shapefile as shp
from matplotlib.ticker import StrMethodFormatter

# Open fire weather districts shapefile



import geopandas
fwd_shapefile_path = '/g/data/ct18/dw4060/viirs/inputfiles/boundaries/FWFDistricts_MapID_WGS84_May2022.shp'
#read shape file into geopandas data frame
df = geopandas.read_file(fwd_shapefile_path)
# we currently care about one district - Mallee
# we need to get the lat.lon of al the points in the polygon (geometry)
fwd_geometry = df[df['DIST_NAME'] == 'Mallee'].geometry.values

import datetime
import glob
import os
file_path = '/g/data/ct18/Historical_Curing/BOM_raw_daily_input_for_VISCA/unzipped/vic/'

year = 2015
t0 = datetime.datetime(year, 7, 1)
#define arrays 
time_array = []
mean_curing_value = []

#8 is a step 
for doy in range(0, 365, 8):
    current_date = t0 + datetime.timedelta(days=doy)
    fns = glob.glob(os.path.join(file_path, "curing-filled_modis-mapvictoria_8day_500m_vic_%s.tif" % current_date.strftime('%Y%m%d')))
    print("fns", fns)
    #cheking file exists
    if len(fns) == 1: 
        #opening the file using xr
        ds = xr.open_dataset(fns[0])
        #taking the cuirng data and clipping I with the geometry of the fwd
        gci_clipped_to_fwd = ds.rio.clip(fwd_geometry, df.crs)
        #print(np.shape(gci_clipped_to_fwd))
        #calculte the mean for the clipped gci
        print(doy, current_date, fns[0], gci_clipped_to_fwd.mean())
        #determine mean curing values for the Mallee fwd and structure them as an array starting from the first idex
        mean_curing_value.append(gci_clipped_to_fwd.mean().to_array().values[0])
        #add current date variable to the time arrazy variable 
        time_array.append(current_date)
        print('mean_curing_value', mean_curing_value)
        print('time_array', time_array)
#plotting out array
plt.plot(time_array,mean_curing_value)
plt.show()
plt.savefig("curing_timeline.png")

#next steps: extract landcover shape values, change axis range, add indicators on the plot at the three phases  