#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: esalathe

New version of code to plot a map with WRF data
This uses pcolormesh which is a lot faster than the contour routine and give a better image
"""


from numpy import (
        linspace,array, log,exp,sin,cos,sqrt, pi,e, 
        zeros, ones, amin,amax, argmax, arange
        )
from matplotlib.pyplot import plot, show

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.cm import get_cmap
import cartopy
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature

"""
This is one option for structuring python code
Define a main() module that does the actual work.
This is defined here, but isn't run until the bottom
of the script file.
"""

def main():
    # open netcdf file
    #ncGridFile = Dataset("grid_d02.nc", "r", format="NETCDF4")
    from netCDF4 import Dataset

    ncFile = Dataset("Netcdf_Files/ccsm4-wrf_1970-2099_T2MAXextr.nc",
                           "r", format="NETCDF4")
        
    # get the data to plot
    T90=ncFile.variables["T2MAX90"][:]
    lats = ncFile.variables["XLAT"][:]
    lons = ncFile.variables["XLONG"][:]
    
    
    # set variable specific options
    plotvar=T90[0,:,:]-273.15
    ## min=round(amin(plotvar)/5)*5+5
    ## max=round(amax(plotvar)/5)*5+5
    Tmin=10
    Tmax=40
    
    # draw plot
    
    plt.figure(figsize=(15,10))
    WRFplot(plotvar,lats,lons, Tmin,Tmax, "CCSM4-WRF 1970-2006", "T90 in °C", "RdYlBu_r")
    plt.show()

"""
After main() define subroutines that are called by main()

Structuring it this way helps to keep the variable namespace
distinct between all the parts of the code
"""
    
# Plotting Routine

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


"""
Once everything is defined, execute it
other script commands can be included here
for example handing errors, inputs, global
variables, etc
"""
    
# Execute
main()

