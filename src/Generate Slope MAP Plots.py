# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 11:39:30 2021

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

#script for generating a single image with three panes of Slope plots
    # for each model pair
#         -temperature wise
# 			-For T90:
# 			1) Slope GCM
# 			2) Slope WRF
# 			3) Difference

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
    

def generate_Single_Slope_Plots():
    t_Data_Types = ["tasmax90"]
    p_Data_Types = ["pr95"]
    wrf_Folder = r"Netcdf_Files" + "\\" +"wrf_Netcdf_Files"
    gcm_Folder = r"Netcdf_Files" + "\\" +"gcm_Netcdf_Files"
    #for each wrf and gcm model pair, generate Slope plot for each of the sampled data types
    for wrf_File in os.listdir(wrf_Folder): 
        wrf_head, wrf_tail = os.path.split(wrf_File) #tail gives the file name and type
        wrf_File = "Netcdf_Files" + "\\" + "wrf_Netcdf_Files" + "\\" + wrf_tail
        for gcm_File in os.listdir(gcm_Folder):
            gcm_head, gcm_tail = os.path.split(gcm_File) #tail gives the file name and type
            wrf_Model = wrf_tail.split("-")[0]
            gcm_File = "Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" +  gcm_tail
            if(gcm_tail.__str__().find(wrf_Model)>=0):
                
            
                if(wrf_tail.__str__().find("T2MAX")>=0 and gcm_tail.__str__().find("tasmax")>=0):
                    for temperature_Type in t_Data_Types:
                        wrf_Data_Type = " "
                        if(temperature_Type == "tasmax90"):
                            wrf_Data_Type = "T2MAX90"
                        elif(temperature_Type == "tasmaxx"):
                            wrf_Data_Type = "T2MAXx"
                        elif(temperature_Type !="tasmax90" and temperature_Type != "tasmaxx"):
                            print("invalid temperature_Type")
                            return 1
                        wrf_Data = Dataset(wrf_File, "r" , format = "NETCDF4")
                        gcm_Data = Dataset(gcm_File, "r" , format = "NETCDF4") 
                        wrf_Slope = wrf_Data.variables["regressionValues_Slope_for_" + wrf_Data_Type][:]
                        gcm_Slope = gcm_Data.variables["Interpolated Slope data based on " + temperature_Type][:]
                        wrf_lon = wrf_Data.variables["XLONG"][:] 
                        wrf_lat = wrf_Data.variables["XLAT"][:]
                        wrf_Slope = wrf_Slope[4:-4,4:-4]
                        wrf_lon = wrf_lon[4:-4,4:-4]
                        wrf_lat = wrf_lat[4:-4,4:-4]
                        gcm_Slope = gcm_Slope[4:-4,4:-4]
                       
                        plt.figure(figsize=(7.5,5))
                        WRFplot(wrf_Slope, wrf_lat, wrf_lon,  amin(wrf_Slope),amax(wrf_Slope), "Slope for "  + wrf_File.split("\\")[2].split("-")[0] + "-wrf" 
                                    + " based on " + wrf_Data_Type ,"Slope", "RdYlBu_r")
                        plt.savefig("WRF Slope Plots"+ "\\" + "Slope for "  + wrf_File.split("\\")[2].split("-")[0] + "-wrf" 
                                    + " based on " + wrf_Data_Type + ".png")
                        plt.figure(figsize=(7.5,5))
                        WRFplot(gcm_Slope, wrf_lat, wrf_lon,  amin(gcm_Slope),amax(gcm_Slope), "Slope for "  + gcm_File.split("\\")[2].split("_")[0] + "-gcm" 
                                    + " based on " + temperature_Type ,"Slope", "RdYlBu_r")
                        plt.savefig("GCM Slope Plots"+ "\\" + "Slope for "  + gcm_File.split("\\")[2].split("_")[0] + "-gcm" 
                                    + " based on " + temperature_Type + ".png")
                        wrf_Data.close()
                        gcm_Data.close()
                        
                        
                if(wrf_tail.__str__().find("PREC")>=0 and gcm_tail.__str__().find("pr_extr")>=0):
                    for precipitation_Type in p_Data_Types:
                        wrf_Data_Type = " "
                        if(precipitation_Type == "pr95"):
                            wrf_Data_Type = "PREC95"
                        elif(precipitation_Type == "prx"):
                            wrf_Data_Type = "PRECx"
                        elif(precipitation_Type != "prx" and precipitation_Type != "pr95"):
                            print("invalid precipitation_Type")
                            return 1
                        wrf_Data = Dataset(wrf_File, "r" , format = "NETCDF4")
                        gcm_Data = Dataset(gcm_File, "r" , format = "NETCDF4") 
                        wrf_Slope = wrf_Data.variables["regressionValues_Slope_for_" + wrf_Data_Type][:]
                        gcm_Slope = gcm_Data.variables["Interpolated Slope data based on " + precipitation_Type][:]
                        wrf_lon = wrf_Data.variables["XLONG"][:] 
                        wrf_lat = wrf_Data.variables["XLAT"][:]
                        wrf_Slope = wrf_Slope[4:-4,4:-4]
                        wrf_lon = wrf_lon[4:-4,4:-4]
                        wrf_lat = wrf_lat[4:-4,4:-4]
                        gcm_Slope = gcm_Slope[4:-4,4:-4]
                        plt.figure(figsize=(7.5,5))
                        WRFplot(wrf_Slope, wrf_lat, wrf_lon,  amin(wrf_Slope),amax(wrf_Slope), "Slope for "  + wrf_File.split("\\")[2].split("-")[0] + "-wrf" 
                                    + " based on " + wrf_Data_Type ,"Slope", "RdYlBu_r")
                        plt.savefig("WRF Slope Plots"+ "\\" + "Slope for "  + wrf_File.split("\\")[2].split("-")[0] + "-wrf" 
                                    + " based on " + wrf_Data_Type + ".png")
                        plt.figure(figsize=(7.5,5))
                        WRFplot(gcm_Slope, wrf_lat, wrf_lon,  amin(gcm_Slope),amax(gcm_Slope), "Slope for "  + gcm_File.split("\\")[2].split("_")[0] + "-gcm" 
                                    + " based on " + precipitation_Type ,"Slope", "RdYlBu_r")
                        plt.savefig("GCM Slope Plots"+ "\\" + "Slope for "  + gcm_File.split("\\")[2].split("_")[0] + "-gcm" 
                                    + " based on " + precipitation_Type + ".png")
                        wrf_Data.close()
                        gcm_Data.close()
                        
#execute
generate_Single_Slope_Plots()