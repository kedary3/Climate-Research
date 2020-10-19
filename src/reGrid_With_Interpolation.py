#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 17:49:02 2020

@author: salathe
"""

####

#to dos
    # make the interpolation x,y,x2,y2 be dynamic data.variables[lat][:]
    # make the interpolation a method call and populate interpolated toe data into gcm files 
    # boundary will be fixed 
    # try to make a toe delta graph (wrf-gcm)
    # find extremes for different variavbles 
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
import pylab as lb
from numpy import *
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

data=Dataset("Netcdf_Files/CCSM4_1970-2099_prextr.nc", "r", format="NETCDF4")
pr95=data.variables["pr95"][:] 
# Here use x as lon and y as lat and Z as data field from GCM
x = linspace(0, 1, 33)
y = linspace(0,1,21)
Z = pr95[0,:,:]

# here X2 and Y2 (upper case) are lon and lat from WRF grid.nc
# note that lat and lon are 2d arrays 
x2 = linspace(0, 1, 123)
y2 = linspace(0, 1, 162)
X2,Y2 = meshgrid(x2,y2) # from grid.nc

f = interp2d(x, y, Z, kind='cubic')

# Since X2 and Y2 are 2-d arrays for irregular grid need to do point by point
#    **THERE MAY BE A BETTER WAY?
Z2=zeros_like(X2)
for i in arange(X2.shape[0]):
    for j in arange(X2.shape[1]):
        Z2[i,j]=f(X2[i,j],Y2[i,j])

# plot coarse and fine grids
X, Y = meshgrid(x, y)

fig, ax = plt.subplots(nrows=1, ncols=2)
ax[0].pcolormesh(X, Y, Z, shading="auto")

ax[1].pcolormesh(X2, Y2, Z2, shading="gouraud")

plt.show()