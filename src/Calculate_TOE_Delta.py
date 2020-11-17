# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 17:31:46 2020

@author: Kedar
"""

#plot a delta for each combination of data type and model
#Delta toe for access and ccsm4 between wrf and gcm for temp and try prec
		#all combinations
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

#plot function
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
    
#define clone to have a suitable file to store delta data
def clone(src_file, trg_file): #function to copy attributes, variables, and dimensions from one NETCDF to another
    src = nc.Dataset(src_file)
    trg = nc.Dataset(trg_file, mode='w')

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions)

        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})

        # Copy the variables values (as 'f4' eventually)
        trg.variables[name][:] = src.variables[name][:]
    # return the file and close the source
    src.close
    return trg
#define helper delta calculator.
def Calculate_TOE_Delta(wrf_File, gcm_File, delta_File, gcm_Data_Type):
    wrf_Data_Type = " "
    if(gcm_Data_Type == "tasmax90"):
        wrf_Data_Type = "T2MAX90"
    elif(gcm_Data_Type == "tasmaxx"):
        wrf_Data_Type = "T2MAXx"
    elif(gcm_Data_Type == "pr95"):
        wrf_Data_Type = "PREC95"
    elif(gcm_Data_Type == "prx"):
        wrf_Data_Type = "PRECx"
    elif(gcm_Data_Type !="tasmax90" and gcm_Data_Type != "tasmaxx" and gcm_Data_Type != "prx" and gcm_Data_Type != "pr95"):
        print("invalid gcm_Data_Type")
        return 1
    wrf_Data = Dataset(wrf_File, "r" , format = "NETCDF4")
    gcm_Data = Dataset(gcm_File, "r" , format = "NETCDF4")
    delta_Data = Dataset(delta_File, "r+", format = "NETCDF4")
    wrf_TOE = wrf_Data.variables["ToE_for_" + wrf_Data_Type][:]
    gcm_TOE = gcm_Data.variables["Interpolated ToE data based on " + gcm_Data_Type][:]
    if(delta_Data.variables.__str__().find("Delta between "  + wrf_File.split("\\")[1].split("-")[0] + "-wrf" +
                                      " and " + gcm_File.split("\\")[1].split("_")[0] + "-gcm" + " based on " + gcm_Data_Type)==-1):
        delta = delta_Data.createVariable("Delta between "  + wrf_File.split("\\")[1].split("-")[0] + "-wrf" +
                                      " and " + gcm_File.split("\\")[1].split("_")[0] + "-gcm" + " based on " + gcm_Data_Type,
                                      "float32", ("south_north", "west_east",))
        delta.units = "years"
    deltas = delta_Data.variables["Delta between "  + wrf_File.split("\\")[1].split("-")[0] + "-wrf" +
                                      " and " + gcm_File.split("\\")[1].split("_")[0] + "-gcm" + " based on " + gcm_Data_Type]
    
    #for each grid point assign a delta based on (diff=wrf-gcm)
    ny = shape(wrf_TOE)[0] # south-north - 123 - wrf / south-north - 21 - gcm
    nx = shape(wrf_TOE)[1] # west-east - 162 - wrf / west-east - 33 - gcm
    for k in range(ny):
        for i in range(nx):
            d = wrf_TOE[[k],[i]]-gcm_TOE[[k],[i]]
            deltas[[k],[i]]=d
    
    wrf_lon=delta_Data.variables["XLONG"][:] 
    wrf_lat=delta_Data.variables["XLAT"][:]
    plt.figure(figsize=(15,10))
    WRFplot(deltas, wrf_lat, wrf_lon,  amin(deltas),amax(deltas), "Delta between "  + wrf_File.split("\\")[2].split("-")[0] + "-wrf" +
                                      " and " + gcm_File.split("\\")[2].split("_")[0] + "-gcm" + " based on " + gcm_Data_Type ,"Difference in TOE in Years", "RdYlBu_r")
    plt.savefig("Delta between "  + wrf_File.split("\\")[2].split("-")[0] + "-wrf" +
                                      " and " + gcm_File.split("\\")[2].split("_")[0] + "-gcm" + " based on " + gcm_Data_Type + ".png")
    wrf_Data.close()
    gcm_Data.close()
    
    delta_Data.close()
    
#execute

#define netcdf to store delta data
if(os.path.exists("Netcdf_Files" + "\\" + "TOE_Deltas_from_gcm_to_wrf.nc") == False):
    delta_File = "Netcdf_Files" + "\\" + "TOE_Deltas_from_gcm_to_wrf.nc"
    ds=clone("Netcdf_Files" + "\\" + "wrf_grid.nc", delta_File)

#for each wrf and gcm file, generate toe data for each data type and get interpolated gcm data
delta_File = "Netcdf_Files" + "\\" + "TOE_Deltas_from_gcm_to_wrf.nc"
# (wrf, gcm) (data_Type)
t_Data_Types = ["tasmax90", "tasmaxx"]
p_Data_Types = ["pr95", "prx"]
wrf_Folder = r"Netcdf_Files" + "\\" +"wrf_Netcdf_Files"
gcm_Folder = r"Netcdf_Files" + "\\" +"gcm_Netcdf_Files"
 
#for each wrf and gcm file, generate toe data for each data type
for wrf_File in os.listdir(wrf_Folder): 
    wrf_head, wrf_tail = os.path.split(wrf_File) #tail gives the file name and type
    for gcm_File in os.listdir(gcm_Folder):
        gcm_head, gcm_tail = os.path.split(gcm_File) #tail gives the file name and type
        if(wrf_tail.__str__().find("T2MAXextr")>=0 and gcm_tail.__str__().find("tasmaxextr")>=0):
            for temperature_Type in t_Data_Types:
                Calculate_TOE_Delta("Netcdf_Files" + "\\" + "wrf_Netcdf_Files" + "\\" + wrf_tail, "Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" +  gcm_tail, delta_File, temperature_Type)  
        if(wrf_tail.__str__().find("PRECextr")>=0 and gcm_tail.__str__().find("prextr")>=0):
            for precipitation_Type in p_Data_Types:
                Calculate_TOE_Delta("Netcdf_Files" + "\\" + "wrf_Netcdf_Files" + "\\" + wrf_tail, "Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" +  gcm_tail, delta_File, precipitation_Type)  
                
