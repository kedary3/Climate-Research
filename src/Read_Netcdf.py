#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: esalathe
"""


from numpy import (
        linspace,array, log,exp,sin,cos,sqrt, pi,e, 
        zeros, ones, amin,amax, argmax, arange, shape,
        floor
    )
import numpy as np

from netCDF4 import Dataset
# https://unidata.github.io/netcdf4-python/netCDF4/index.html

# open netcdf file
ncFile = Dataset("ccsm4_1970-2006_TX90th.nc",
                       "r", format="NETCDF4")

# Get the data from the netcdf file
#
# Note that using these two forms gives different results:
# t90 = ncFile.variables["T90"]
# t90 = ncFile.variables["T90"][:]
#
# The first form creates an object with all the metadata that goes with the variable t90
# The [:] form forces strips the metadata and leaves only the bare numpy array.
# Not sure this is a documented feature, but everyone knows it...

t90 = ncFile.variables["T90"][:] - 273.15 # convert to deg-C
lat = ncFile.variables["XLAT"][:]
lon = ncFile.variables["XLONG"][:]

# get dimensions
ntime = shape(t90)[0] # year index
ny = shape(t90)[1] # south-north
nx = shape(t90)[0] # west-east

# make an array with the year values
year0 = 1970
year = arange(ntime)+year0

# find gridpoint closest to an arbitrary lat-lon point
lat_point = 45.
lon_point = -122.

xx = sqrt( (lat - lat_point)**2   +   (lon - lon_point)**2)
idx = np.where( xx == xx.min())
jy=idx[0]
ix=idx[1]

print("Latitude:  ", lat[jy,ix])
print("Longitude: ", lon[jy,ix])

# plot the data at that point

from pylab import *
figure()
plot(year,t90[:,jy,ix])
title("Annual 90th Percentile Tmax")
ylabel("Tmax (deg-C)")
xlabel("Year")
savefig("t90.png", bbox_inches='tight')
close()

        