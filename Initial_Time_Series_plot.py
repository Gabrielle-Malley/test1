#%%
from netCDF4 import Dataset as NetCDFFile
#mapping libraries
from matplotlib.cm import get_cmap
import cartopy.crs as ccrs
import cartopy
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import AxesGrid
#for reading the shape file

from rasterio import features
from affine import Affine
#read subroutines from another file
from shape_routines import add_shape_coord_from_data_array
from shape_routines import transform_from_latlon
from shape_routines import rasterize

#read in libraries
#standard libraries
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
#read libraries

#Module Name plot_netcdf_map_example
#Date of creation: 24/3/2022
#Name of creator of module : CP
#History of modification :
#Summary of what the module does
# this routine reads in a netcdf file
#typical usage
#typical call
#read_modis_curing(file_path= '/scratch/er8/ljm548/curing_modis-mapvictoria_8day_500m_aust_20210605-20210612.nc',varin='curing',unit_test=1):



# things you might find helpful
# use ncdump -h (on command line)  to look at the header of the netcdf file and to work out what variables are in the file and what the dimensions of the variable are.
# Is it a reglar latitude/longitude grid?
# what are the order or lat and lon in the 2d array
#good resource
#http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html




def read_modis_curing(file_path= '/scratch/er8/ljm548/curing_modis-mapvictoria_8day_500m_aust_20210605-20210612.nc',varin='curing',unit_test=1):

    #check file exists

    if os.path.exists(file_path):
        print('files exists',file_path)
    else:
        print('file does not exist',file_path)


    #read in the data nto a variable called nc
    nc = NetCDFFile(file_path)

    # what is in the file ?writes out variables in the netcdf file
    if unit_test>0:
        print(nc.variables.keys()) 


    #read in the location data # note than in many files it if often 'lat' rather than 'latitude'
    lat = nc.variables['latitude'][:]
    lon = nc.variables['longitude'][:]


    #read in the variable you want to plot
    variable = nc.variables[varin][:,:]
    units=nc.variables[varin].units
    if unit_test>0:
        print('units',units)
        print('check the shape of the variable',np.shape(variable))

    return variable,lon,lat,units,varin




#Module Name plot_netcdf_map_example
#Date of creation: 24/3/2022
#Name of creator of module : CP
#History of modification :
#Summary of what the module does
# this routine plots a map from input 2D varable and 1D lat lon
#typical usage
#read_modis_curing(file_path= '/scratch/er8/ljm548/curing_modis-mapvictoria_8day_500m_aust_20210605-20210612.nc',varin='curing',unit_test=1):


def plot_curing(variable,lon,lat,units,varin,plotdir='/g/data/er8/global_challenge/code/plots/',region=22,unit_test=1,avlat=None,avlon=None):


        #get central lat lon of map
    if avlat==None:
        mlat=np.mean(lat)
        mlon=np.mean(lon)
    else:
        mlat=avlat
        mlon=avlon

    #sets a region around the central point to plot
    minlat=mlat-region
    maxlat=mlat+region

    minlon=mlon-region
    maxlon=mlon+region

    #setting the projection to plot
    #https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html
    ax = plt.axes(projection = ccrs.PlateCarree())

    #set the region to plot
    ax.set_extent([minlon, maxlon, minlat, maxlat], crs = ccrs.PlateCarree())
    ax.coastlines()

    #setting min/max values of the plot
    minv=np.min(variable)
    maxv=np.max(variable)


    #mesh implies it only need idimensional lat/lon variables but a 2D variable parameter.
    #https://matplotlib.org/3.5.1/tutorials/colors/colormaps.html
    #lots of different colour masks
    #use '_r' to reverse the color scale
    
    im = ax.pcolormesh(lon, lat, np.squeeze(variable), vmin = minv, vmax = maxv, cmap = 'RdYlGn_r')

    #plot the color bar
    cb = plt.colorbar(im, pad = 0.1, orientation = 'horizontal', shrink = 0.6)
    cb.set_label(units)

    # make the map look nice Plot gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels = True, color = 'gray', linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Create a title and show the data
    title = varin.upper()

    plt.title(title)
    #save map to a file
    plt.savefig(plotdir+title +'.png')
    #plot the map o the screen
    if unit_test>0:
        plt.show()
    plt.close()


###now call the above routines

#read the data        
variable,lon,lat,units,varin=read_modis_curing(file_path= '/scratch/er8/ljm548/curing_modis-mapvictoria_8day_500m_aust_20210605-20210612.nc',varin='curing',unit_test=1)



#plot
plot_curing(variable,lon,lat,units,varin,plotdir='/g/data/er8/global_challenge/code/plots/',unit_test=1,region=22)