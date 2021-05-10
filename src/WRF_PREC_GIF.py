# -*- coding: utf-8 -*-
"""
Created on Mon May 10 04:30:41 2021

@author: Kedar
"""

from numpy import (
        linspace,array, log,exp,sin,cos,sqrt, pi,e, 
        zeros, ones, amin,amax, argmax, arange, shape
        )
from matplotlib.pyplot import plot, show

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.cm import get_cmap
import cartopy
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
import pylab as lb
from netCDF4 import Dataset
import netCDF4 as nc
import pylab as lb
from numpy import *
from scipy.interpolate import interp2d
import os
import numpy.ma as ma

import os
import numpy as np
import matplotlib.pyplot as plt
import imageio
#precipitation gif 

def WRFplot(plotvar,lats,lons,  vmin,vmax, title,varname, ColMap):
    # Set the cartopy mapping object for the WRF domain
    #  (Taken from wrf getcartopy and cartopy_xlim)
    cart_proj = cartopy.crs.LambertConformal(
            central_longitude=-121.0, 
            central_latitude=45.665584564208984, 
            false_easting=0.0, 
            false_northing=0.0, 
            secant_latitudes=None, 
            standard_parallels=[30.,60.], 
            globe=None, 
            cutoff=-30)
    ax = plt.axes(projection=cart_proj)
    ax.set_xlim([-875806.9669240027, 1056192.549175313])
    ax.set_ylim([-733768.6404772081, 730230.3670079684])
    
    # Download and add the states and coastlines
    states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_1_states_provinces_lines")
    ax.add_feature(states, linewidth=.5, edgecolor="black")
    borders = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_0_boundary_lines_land")
    ax.add_feature(borders, linewidth=.75, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)
    
    # Color in the data on the map with smoothing
    plt.pcolormesh(lons,lats,
                   plotvar,vmin=vmin,vmax=vmax,
                   transform=crs.PlateCarree(),
                   shading='gouraud',
                   cmap=get_cmap(ColMap)
                   )
    
    # Add a color bar
    cbar=plt.colorbar(ax=ax, shrink=.6)#, orientation='horizontal')
    cbar.set_label(varname)
    
    # Add gridlines
    #ax.gridlines(color="black", linestyle="dotted")
    
    # Add a title
    plt.title(title)   
    

    
def create_WRF_PREC_Gif(file):
    ncFile = Dataset(file,"r", format="NETCDF4")
    PREC=ncFile.variables["PREC95"][:]
    #get time axis and grid values
    year=shape(PREC)[0]
    years=130
    ny = shape(PREC)[1] # south-north - 123
    nx = shape(PREC)[2] # west-east - 162
    wrf_lon=ncFile.variables["XLONG"][:] 
    wrf_lat=ncFile.variables["XLAT"][:]
    PREC = PREC[:,4:-4,4:-4]
    wrf_lon = wrf_lon[4:-4,4:-4]
    wrf_lat = wrf_lat[4:-4,4:-4]
    
    text_kwargs = dict(ha='left', va='bottom', fontsize=28, color='C1')
    left, width = .25, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height
    

    filenames = []
    for i in np.arange(years-1):
        # plot the year's precipitation
        plt.figure(figsize=(15,10))
        WRFplot(PREC[i][:][:], wrf_lat, wrf_lon,  amin(PREC),amax(PREC), "95th Percentile Precipitation Based on Access1.3 Over Time", "PREC95 (Kg per m^2)", "RdYlBu_r")
        plt.text(left, bottom, "Year: "+str(i),
        horizontalalignment='right',
        verticalalignment='bottom', **text_kwargs)
        
        plt.savefig("year" + str(i) + ".png")
        
        
        # create file name and append it to a list
        filename = "year" + str(i) + ".png"
        filenames.append(filename)
        
        # save frame
        plt.savefig(filename)
        plt.close()
    # build gif
    with imageio.get_writer('mygif.gif', mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)
            
    # Remove files
    for filename in set(filenames):
        os.remove(filename)
        
    ncFile.close()
create_WRF_PREC_Gif("Netcdf_Files\\wrf_Netcdf_Files\\access1.3-wrf_1970-2099_PREC_extr.nc")