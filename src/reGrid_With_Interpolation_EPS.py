#!/opt/local/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 17:49:02 2020

@author: salathe
"""

#to dos
    # make the interpolation x,y,x2,y2 be dynamic data.variables[lat][:]
    # make the interpolation a method call and populate interpolated toe data into gcm files
    # boundary will be fixed
    # try to make a toe delta graph (wrf-gcm)
    # find extremes for different variables 

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



data=Dataset("Netcdf_Files/CCSM4_1970-2099_prextr.nc", "r", format="NETCDF4")
pr95=data.variables["pr95"][:] * 24*60*60 # mm/day
lon=data.variables["lon"][:] - 360
lat=data.variables["lat"][:] 

# Here use x as lon and y as lat and Z as data field from GCM
x = lon
y = lat
Z = pr95[0,:,:]

# here X2 and Y2 (upper case) are lon and lat from WRF grid.nc
# note that lat and lon are 2d arrays 
#x2 = linspace(0, 1, 123)
#y2 = linspace(0, 1, 162)
#X2,Y2 = meshgrid(x2,y2) # from grid.nc

wdata=Dataset("Netcdf_Files/wrf_grid.nc", "r", format="NETCDF4")
wrf_lon=wdata.variables["XLONG"][:] 
wrf_lat=wdata.variables["XLAT"][:] 


wrf_interp = interp2d(lon, lat, pr95[0,:,:], kind='cubic')

# Since X2 and Y2 are 2-d arrays for irregular grid need to do point by point
#    **THERE MAY BE A BETTER WAY?
wrf_pr=zeros_like(wrf_lon)
for i in arange(wrf_lon.shape[0]):
    for j in arange(wrf_lon.shape[1]):
        wrf_pr[i,j]=wrf_interp(wrf_lon[i,j],wrf_lat[i,j])

# plot coarse and fine grids
LON, LAT = meshgrid(lon,lat)

fig, ax = plt.subplots(nrows=1, ncols=2)
ax[0].pcolormesh(LON,LAT,[0,:,:], shading="auto")

ax[1].pcolormesh(wrf_lon,wrf_lat,wrf_pr, shading="auto")


plt.figure(figsize=(15,10))
WRFplot(wrf_pr,wrf_lat,wrf_lon, amin(wrf_pr),amax(wrf_pr), "CCSM4 1970-2006", "PR max in mm/day", "RdYlBu_r")
plt.show()


plt.show()



