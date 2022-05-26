#this script from Dani Martin CFA
#this code reads in either a cfa curing file in tif format
#or a netcdf curing data
#It extracts la and lon information form the tif file using rasterio
#it can also read in a shape files, e.g districts states etc to mask a region
# takes a little while to run
# import libraries
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

fwd_shapefile_path = '/g/data/ct18/dw4060/viirs/inputfiles/boundaries/BOM_FireForecastDistricts_Australia_WGS1984.shp'
#read in the shape file
fwd_shapefile = shp.Reader(fwd_shapefile_path)

# select jurisdictions
cfa_flag=1 # set this to one when reading a tif file to 0 when reading a netcdf file
jurisdictions = ["vic"]
# eg. jurisdictions = ["vic", "tas", "nsw", "qld", "sa", "act"]

for jurisdiction in jurisdictions:

    #Define id numbers of districts
    
    if jurisdiction == 'vic':
        i_list = [27,32,35,36,37,38,39,40,41]

    # iterate through dates
    if cfa_flag>0:
        dates = ["20190223", "20190224"]
    else:
        dates = ["20210605-20210612"]
        

    for date in dates:
        
        # Define filenames and directories
        if cfa_flag==1:
            modis_folder = '/g/data/ct18/Historical_Curing/BOM_raw_daily_input_for_VISCA/unzipped/'+jurisdiction
            extension='.tif'

            
        else:
            modis_folder = 'netcdf:/scratch/er8/ljm548/curing_modis-mapvictoria_8day_500m_aust_20210605-20210612.nc:curing'
            extension='.nc'
            
            

        # define modis filename (filenames vary between jurisdictions)
        
        if jurisdiction == 'act':
            modis_filename = 'curing-filled_modis-mapvictoria_8day_500m_'+'vic'+'_'+date+extension
        elif jurisdiction == 'vic':
            modis_filename = 'curing-filled_modis-mapvictoria_8day_500m_'+jurisdiction+'_'+date+'.legend'+extension
        else:
            modis_filename = 'curing-filled_modis-mapvictoria_8day_500m_'+jurisdiction+'_'+date+extension
        if cfa_flag>0:
            modis_filepath = os.path.join(modis_folder, modis_filename)
        else:
            modis_folder='/scratch/er8/ljm548/'

            modis_filename='curing_modis-mapvictoria_8day_500m_aust_'+date+'.nc'
            
            modis_filepath = os.path.join(modis_folder, modis_filename)
            
        for modis_file in glob.glob(modis_filepath):
        
            # Check if modis files exist
            
            if os.path.isfile(modis_file):
                    
                # plot fire weather district shapefile on map
                
                for i in i_list:
                    shape = fwd_shapefile.shape(i)
                    shaperecord = fwd_shapefile.shapeRecord(i)
                    shaperecord_shape = shaperecord.shape
                    points = np.array(shape.points)       
                    intervals = list(shape.parts) + [len(shape.points)]
                    ax = plt.gca()
                    ax.set_aspect(1)
                    for (i, j) in zip(intervals[:-1], intervals[1:]):
                        ax.plot(*zip(*points[i:j]), color='black')

                    # define curing color ramp
                    
                    curingcmap = ListedColormap(["darkgreen", "forestgreen", "limegreen", "greenyellow", "yellow", "gold", "orange", "darkorange", "orangered", "red"])

                    # plot modis curing data on map
                    if cfa_flag>0:
                        modisraster = rio.open_rasterio(modis_file)
                        #type raster <class 'xarray.core.dataarray.DataArray'>

                    else:
                        print('netcdf:'+modis_file+':curing')
                        modis_filex= 'netcdf:'+modis_file+':curing'
                        print('modis_filex',modis_filex)
                        import rioxarray
                        #Leon modisraster = rasterio.open(modis_filex)
                        #type raster <class 'rasterio.io.DatasetReader'>
                        modisraster = rioxarray.open_rasterio(modis_filex)
                        #type raster <class 'xarray.core.dataarray.DataArray'>
                    print('type raster',type(modisraster))
                        

                    modisarray = xr.where(modisraster > 100.00, np.nan, modisraster)
                    #sys.exit()
                    modisarray = xr.where(modisarray < 0.00, np.nan, modisarray)        
                    modisarray.plot(cmap=curingcmap)
                    plt.title(shaperecord.record.DIST_NAME+": MODIS ("+ date+")")
                    plt.xlabel(None)
                    plt.ylabel(None)
                    plt.locator_params(axis='x', nbins=3)
                    plt.locator_params(axis='y', nbins=3)
                    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
                    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
                    print('coords',shape.bbox[0], shape.bbox[2], shape.bbox[1], shape.bbox[3])

                    #extract out single specific lat lon array
                    
                    xs = np.array([shape.bbox[0], shape.bbox[2]])
                    ys = np.array([shape.bbox[1], shape.bbox[3]])
                    if cfa_flag>0:
                    #https://stackoverflow.com/questions/63004971/find-latitude-longitude-coordinates-of-every-pixel-in-a-geotiff-image
                    
                        with rasterio.open(modis_file) as src:
                            rows, cols = rasterio.transform.rowcol(src.transform, xs, ys)
                            print(src.profile)
                            #extract out lon and lat arrays
                            r=src.read(1)
                            print(r.shape)
                            width=len(r[:,0])
                            height=len(r[0,:])

                            print(width,height)
                            latarray=np.zeros([width,height])
                            lonarray=np.zeros([width,height])
                            for iwidth in range(width):
                                for iheight in range(height):

                                    pixels2coords = src.xy(iwidth,iheight)  #input px, py

                                    lonarray[iwidth,iheight]=pixels2coords[0]
                                    latarray[iwidth,iheight]=pixels2coords[1]

                            print('lonarray',lonarray)
                            print(np.shape(lonarray))
                            print(np.shape(modisarray))
                            print(modisarray)
                            print(np.shape(cols))

                
                    plt.axis([shape.bbox[0]-0.05, shape.bbox[2]+0.05, shape.bbox[1]-0.05, shape.bbox[3]+0.05])
                    
                    # save maps (and amend 'Illawarra/Shoalhaven' district name)    
                                   
                    if shaperecord.record.DIST_NAME == 'Illawarra/Shoalhaven':
                        plt.savefig(r'/home/565/gm6737/Trial_Code_Files'+jurisdiction.strip()+'_'+date+'_'+'IllawarraShoalhaven'+'.png')
                    elif shaperecord.record.DIST_NAME != 'Illawarra/Shoalhaven':
#                        plt.savefig(r'/g/data/er8/users/cp2786//global_challenge/outputfiles/examples/Map_MODIS_Curing_'+jurisdiction.strip()+'_'+date+'_'+shaperecord.record.DIST_NAME.strip()+'.png')


                        filenameout='/home/565/gm6737/Trial_Code_Files'+jurisdiction.strip()+'_'+date+'_'+shaperecord.record.DIST_NAME.strip()
                        filenameout=filenameout.strip()
                        plt.savefig(filenameout+'.png')

                    else:
                        filenameout='/home/565/gm6737/Trial_Code_Files'+jurisdiction.strip()+'_'+date+'_'+shaperecord.record.DIST_NAME.strip()
                        filenameout=filenameout.strip()
                        plt.savefig(filenameout+'.png')
                     
                    plt.close()                           


