#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: esalathe
"""

# This reads in a netcdf file on the WRF model domain and plots tjhe data on a map
# the cartography uses the python cartopy package
# https://scitools.org.uk/cartopy/docs/latest/

from numpy import (
        linspace,array, log,exp,sin,cos,sqrt, pi,e, 
        zeros, ones, amin,amax, argmax, arange
        )
from matplotlib.pyplot import plot, show

import numpy as np
import os

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature


# Plotting Routine

def WRFplot(plotvar,plotvarname):
    # Set the GeoAxes to the projection used by WRF
    ax = plt.axes(projection=cart_proj)
    
# Download and add the states and coastlines
    states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_1_states_provinces_lines")
    ax.add_feature(states, linewidth=.5, edgecolor="black")
    borders = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_0_boundary_lines_land")
    ax.add_feature(borders, linewidth=1, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)
    
    # add the filled contour plot
    plt.contourf(lons, lats, plotvar, plt_range,
         transform=crs.PlateCarree(),
         extend=extend,
         cmap=get_cmap(ColMap))

    
    # Add a color bar
    cbar=plt.colorbar(ax=ax, shrink=.6)#, orientation='horizontal')
    cbar.set_label('deg-C')
    
    # add title    
    plt.title(plotvarname)
        

# open netcdf file

ncFile = Dataset("ccsm4_1970-2006_TX90th.nc",
                       "r", format="NETCDF4")

# Set the cartopy mapping object (can extract from netcdf)
# https://scitools.org.uk/cartopy/docs/latest/crs/projections.html
cart_proj = cartopy.crs.LambertConformal(
        central_longitude=-121.0, 
        central_latitude=45.665584564208984, 
        false_easting=0.0, 
        false_northing=0.0, 
        secant_latitudes=None, 
        standard_parallels=[30.,60.], 
        globe=None, 
        cutoff=-30)

# get the data to plot
T90=ncFile.variables["T90"][:]
lats = ncFile.variables["XLAT"][:]
lons = ncFile.variables["XLONG"][:]

# set variable to plot
plotvar=np.mean(T90[:,:,:], axis=0) - 273.15

# set parameters for the plot
min=10
max=40
plt_range=arange(min,max,1)
extend='both'
# https://matplotlib.org/tutorials/colors/colormaps.html
ColMap="coolwarm"


# draw plot

fig = plt.figure(figsize=(15,10))

WRFplot(plotvar,"T90")

plt.show()