# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 20:15:28 2021

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


#use np.nanmean

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
    
def mask(ar): #mask values above the cuttoff
    for i in range(shape(ar)[0]):
        for k in range(shape(ar)[1]):
            if ar[i][k]>2299:
                ar[i][k] = np.nan
    return ar    
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


def Calculate_gcm_Average_TOE_of_Model_Types(avg_File):
    
    models = ["access1.3","bcc","canesm2","ccsm4","csiro","fgoals","gfdl","giss","miroc5","mri","noresm1"]
    data_Types_File_Names = ["tasmax","pr"]
    
    data_Arr = np.empty((11,123,162))
    
    avg = Dataset(avg_File, "r+", format = "NETCDF4")
    
    for data_Type_File_Name in data_Types_File_Names:
        for model in models:
            model_Index = models.index(model)
            
            for gcm_File in os.listdir(gcm_Folder):
                gcm_head, gcm_tail = os.path.split(gcm_File)
                if(gcm_tail.__str__().find(model)>=0 and gcm_tail.__str__().find(data_Type_File_Name)>=0):
                    if(data_Type_File_Name == "pr"):
                        data_Type = "pr95"
                    if(data_Type_File_Name == "tasmax"):
                        data_Type = "tasmax90"
                    gcm_File= "Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" +  gcm_tail
                    gcm_Data = Dataset(gcm_File, "r" , format = "NETCDF4")
                    gcm_TOE = gcm_Data.variables["Interpolated ToE data based on " + data_Type][:]
                    # gcm_TOE = mask(gcm_TOE)
                    
                    data_Arr[model_Index][:][:] = gcm_TOE
                    gcm_Data.close()
                    
        if(avg.variables.__str__().find("Average TOE of gcm models based on " + data_Type)==-1):
            gcm_average_TOE = avg.createVariable("Average TOE of gcm models based on " + data_Type,
                          "float32", ("south_north", "west_east",))
        if(avg.variables.__str__().find("SDV of Average TOE of gcm models based on " + data_Type)==-1):
            gcm_average_TOE_SDV = avg.createVariable("SDV of Average TOE of gcm models based on " + data_Type,
                                  "float32", ("south_north", "west_east",))
        gcm_average_TOE = avg.variables["Average TOE of gcm models based on " + data_Type]
        gcm_average_TOE_SDV = avg.variables["SDV of Average TOE of gcm models based on " + data_Type]
        gcm_TOE_Avg = np.nanmean(data_Arr,axis=0)
        gcm_TOE_SDV = np.std(data_Arr,axis=0)
        ny = shape(gcm_TOE_Avg)[0] # south-north - 123 - wrf / south-north - 21 - gcm
        nx = shape(gcm_TOE_Avg)[1] # west-east - 162 - wrf / west-east - 33 - gcm
        for k in range(ny):
            for i in range(nx):
                
                gcm_average_TOE[[k],[i]] = gcm_TOE_Avg[[k],[i]]
                gcm_average_TOE_SDV[[k],[i]] = gcm_TOE_SDV[[k],[i]]
        
        
    avg.close()                

def Calculate_wrf_Average_TOE_of_Model_Types(avg_File):
    
    models = ["access1.3","bcc","canesm2","ccsm4","csiro","fgoals","gfdl","giss","miroc5","mri","noresm1"]
    data_Types_File_Names = ["T2MAX","PREC"]
    
    data_Arr = np.empty((11,123,162))
    
    avg = Dataset(avg_File, "r+", format = "NETCDF4")
    
    for data_Type_File_Name in data_Types_File_Names:
        for model in models:
            model_Index = models.index(model)
            
            for wrf_File in os.listdir(wrf_Folder):
                wrf_head, wrf_tail = os.path.split(wrf_File)
                if(wrf_tail.__str__().find(model)>=0 and wrf_tail.__str__().find(data_Type_File_Name)>=0):
                    if(data_Type_File_Name == "PREC"):
                        data_Type = "PREC95"
                    if(data_Type_File_Name == "T2MAX"):
                        data_Type = "T2MAX90"
                    
                    wrf_File= "Netcdf_Files" + "\\" + "wrf_Netcdf_Files" + "\\" +  wrf_tail
                    wrf_Data = Dataset(wrf_File, "r" , format = "NETCDF4")
                    wrf_TOE = wrf_Data.variables["ToE_for_" + data_Type][:]
                    # wrf_TOE = mask(wrf_TOE)
                    
                    
                    data_Arr[model_Index][:][:] = wrf_TOE
                    wrf_Data.close()
        
        if(avg.variables.__str__().find("Average TOE of wrf models based on " + data_Type)==-1):
            wrf_average_TOE = avg.createVariable("Average TOE of wrf models based on " + data_Type,
                                      "float32", ("south_north", "west_east",))
        if(avg.variables.__str__().find("SDV of Average TOE of wrf models based on " + data_Type)==-1):
            wrf_average_TOE_SDV = avg.createVariable("SDV of Average TOE of wrf models based on " + data_Type,
                                  "float32", ("south_north", "west_east",))
        wrf_average_TOE = avg.variables["Average TOE of wrf models based on " + data_Type]
        wrf_average_TOE_SDV = avg.variables["SDV of Average TOE of wrf models based on " + data_Type]            
        
        wrf_TOE_Avg = np.nanmean(data_Arr,axis=0)
        wrf_TOE_SDV = np.nanstd(data_Arr,axis=0)
        ny = shape(wrf_TOE_Avg)[0] # south-north - 123 - wrf / south-north - 21 - gcm
        nx = shape(wrf_TOE_Avg)[1] # west-east - 162 - wrf / west-east - 33 - gcm
        for k in range(ny):
            for i in range(nx):
                
                wrf_average_TOE[[k],[i]] = wrf_TOE_Avg[[k],[i]]
                wrf_average_TOE_SDV[[k],[i]] = wrf_TOE_SDV[[k],[i]]
    
    avg.close()                

def Calculate_Difference_Between_Averages(avg_File):
    t_Data_Types = ["tasmax90"]
    p_Data_Types = ["pr95"]
    for temperature_Type in t_Data_Types:
        gcm_Data_Type = temperature_Type
        if(gcm_Data_Type == "tasmax90"):
            wrf_Data_Type = "T2MAX90"
        elif(gcm_Data_Type == "tasmaxx"):
            wrf_Data_Type = "T2MAXx"
        elif(gcm_Data_Type == "pr95"):
            wrf_Data_Type = "PREC95"
        elif(gcm_Data_Type == "prx"):
            wrf_Data_Type = "PRECx"
        elif(gcm_Data_Type !="tasmax90" and gcm_Data_Type != "tasmaxx" and 
             gcm_Data_Type != "prx" and gcm_Data_Type != "pr95"):
            print("invalid gcm_Data_Type")
            return 1
        avg_Data = Dataset(avg_File, "r+" , format = "NETCDF4")
        wrf_average_TOE = avg_Data.variables["Average TOE of wrf models based on " + wrf_Data_Type][:]
        gcm_average_TOE = avg_Data.variables["Average TOE of gcm models based on " + gcm_Data_Type][:]
        if(avg_Data.variables.__str__().find("Difference Between Average WRF and GCM TOE Based On " + gcm_Data_Type)==-1):
                diff_Avg = avg_Data.createVariable("Difference Between Average WRF and GCM TOE Based On " + gcm_Data_Type,
                                          "float32", ("south_north", "west_east",))
        diff_Avgs = avg_Data.variables["Difference Between Average WRF and GCM TOE Based On " + gcm_Data_Type]
        ny = shape(wrf_average_TOE)[0] # south-north - 123 - wrf / south-north - 21 - gcm
        nx = shape(wrf_average_TOE)[1] # west-east - 162 - wrf / west-east - 33 - gcm
        for k in range(ny):
            for i in range(nx):
                
                d = wrf_average_TOE[[k],[i]]-gcm_average_TOE[[k],[i]]
                diff_Avgs[[k],[i]] = d
        avg_Data.close()
        
    for precipitation_Type in p_Data_Types:
        gcm_Data_Type = precipitation_Type
        if(precipitation_Type == "pr95"):
            wrf_Data_Type = "PREC95"
        elif(precipitation_Type == "prx"):
            wrf_Data_Type = "PRECx"
        elif(precipitation_Type != "prx" and precipitation_Type != "pr95"):
            print("invalid precipitation_Type")
            return 1
        avg_Data = Dataset(avg_File, "r+" , format = "NETCDF4")
        wrf_average_TOE = avg_Data.variables["Average TOE of wrf models based on " + wrf_Data_Type][:]
        gcm_average_TOE = avg_Data.variables["Average TOE of gcm models based on " + gcm_Data_Type][:]
        if(avg_Data.variables.__str__().find("Difference Between Average WRF and GCM TOE Based On " + gcm_Data_Type)==-1):
                diff_Avg = avg_Data.createVariable("Difference Between Average WRF and GCM TOE Based On " + gcm_Data_Type,
                                          "float32", ("south_north", "west_east",))
        diff_Avgs = avg_Data.variables["Difference Between Average WRF and GCM TOE Based On " + gcm_Data_Type]
        ny = shape(wrf_average_TOE)[0] # south-north - 123 - wrf / south-north - 21 - gcm
        nx = shape(wrf_average_TOE)[1] # west-east - 162 - wrf / west-east - 33 - gcm
        for k in range(ny):
            for i in range(nx):
                
                d = wrf_average_TOE[[k],[i]]-gcm_average_TOE[[k],[i]]
                diff_Avgs[[k],[i]] = d
        avg_Data.close()

def plot_Averages(avg_File):
    t_Data_Types = ["tasmax90"]
    p_Data_Types = ["pr95"]
    average_Plots_Folder = r"Average TOE of Models Plots"
    for precipitation_Type in p_Data_Types:
        gcm_Data_Type = precipitation_Type
        if(precipitation_Type == "pr95"):
            wrf_Data_Type = "PREC95"
        elif(precipitation_Type == "prx"):
            wrf_Data_Type = "PRECx"
        elif(precipitation_Type != "prx" and precipitation_Type != "pr95"):
            print("invalid precipitation_Type")
            return 1
        avg_Data = Dataset(avg_File, "r" , format = "NETCDF4")
        
        
        avg_lon = avg_Data.variables["XLONG"][:] 
        avg_lat = avg_Data.variables["XLAT"][:]
        avg_WRF_TOE = avg_Data.variables["Average TOE of wrf models based on " + wrf_Data_Type][:]
        avg_GCM_TOE = avg_Data.variables["Average TOE of gcm models based on " + gcm_Data_Type][:]
        avg_Diff = avg_Data.variables["Difference Between Average WRF and GCM TOE Based On " + gcm_Data_Type]
        avg_WRF_TOE_SDV = avg_Data.variables["SDV of Average TOE of wrf models based on " + wrf_Data_Type][:]
        avg_GCM_TOE_SDV = avg_Data.variables["SDV of Average TOE of gcm models based on " + gcm_Data_Type][:]
        avg_WRF_TOE = avg_WRF_TOE[4:-4,4:-4]
        avg_GCM_TOE = avg_GCM_TOE[4:-4,4:-4]
        avg_WRF_TOE_SDV = avg_WRF_TOE_SDV[4:-4,4:-4]
        avg_GCM_TOE_SDV = avg_GCM_TOE_SDV[4:-4,4:-4]
        avg_Diff = avg_Diff[4:-4,4:-4]
        avg_lon = avg_lon[4:-4,4:-4]
        avg_lat = avg_lat[4:-4,4:-4]
        plt.figure(figsize=(7.5,5))
        
        WRFplot(avg_WRF_TOE, avg_lat, avg_lon,  2010,2300,
                "Average TOE of wrf models based on " + wrf_Data_Type,"Time Of Emergence (years)", "RdYlBu_r")
        plt.savefig(average_Plots_Folder + "//" + "Average TOE of wrf models based on " + wrf_Data_Type + ".png")
        plt.figure(figsize=(7.5,5))
        WRFplot(avg_GCM_TOE, avg_lat, avg_lon,  2010,2300,
                "Average TOE of gcm models based on " + wrf_Data_Type,"Time Of Emergence (years)", "RdYlBu_r")
        plt.savefig(average_Plots_Folder + "//" + "Average TOE of gcm models based on " + gcm_Data_Type + ".png")
        
        plt.figure(figsize=(7.5,5))
        WRFplot(avg_Diff, avg_lat, avg_lon,  -1*amax(avg_Diff), amax(avg_Diff),
                "Difference Between Average WRF and GCM TOE Based On " + gcm_Data_Type, "Time Of Emergence (years)", "RdYlBu_r")
        plt.savefig(average_Plots_Folder + "//" + "Difference Between Average WRF and GCM TOE Based On " + gcm_Data_Type+ ".png")
        
        plt.figure(figsize=(7.5,5))
        
        WRFplot(avg_WRF_TOE_SDV, avg_lat, avg_lon,  -1*amax(avg_WRF_TOE_SDV),amax(avg_WRF_TOE_SDV),
                "SDV of Average TOE of wrf models based on " + wrf_Data_Type,"Standard Deviation", "RdYlBu_r")
        plt.savefig(average_Plots_Folder + "//" + "SDV of Average TOE of wrf models based on " + wrf_Data_Type + ".png")
        
        plt.figure(figsize=(7.5,5))
        
        WRFplot(avg_GCM_TOE_SDV, avg_lat, avg_lon,  -1*amax(avg_GCM_TOE_SDV),amax(avg_GCM_TOE_SDV),
                "SDV of Average TOE of gcm models based on " + gcm_Data_Type,"Standard Deviation", "RdYlBu_r")
        plt.savefig(average_Plots_Folder + "//" + "SDV of Average TOE of gcm models based on " + gcm_Data_Type + ".png")
        avg_Data.close()
    for temperature_Type in t_Data_Types:
        gcm_Data_Type = temperature_Type
        if(gcm_Data_Type == "tasmax90"):
            wrf_Data_Type = "T2MAX90"
        elif(gcm_Data_Type == "tasmaxx"):
            wrf_Data_Type = "T2MAXx"
        elif(gcm_Data_Type == "pr95"):
            wrf_Data_Type = "PREC95"
        elif(gcm_Data_Type == "prx"):
            wrf_Data_Type = "PRECx"
        elif(gcm_Data_Type !="tasmax90" and gcm_Data_Type != "tasmaxx" and 
             gcm_Data_Type != "prx" and gcm_Data_Type != "pr95"):
            print("invalid gcm_Data_Type")
            return 1
        avg_Data = Dataset(avg_File, "r" , format = "NETCDF4")
        
        
        avg_lon = avg_Data.variables["XLONG"][:] 
        avg_lat = avg_Data.variables["XLAT"][:]
        avg_WRF_TOE = avg_Data.variables["Average TOE of wrf models based on " + wrf_Data_Type][:]
        avg_GCM_TOE = avg_Data.variables["Average TOE of gcm models based on " + gcm_Data_Type][:]
        avg_Diff = avg_Data.variables["Difference Between Average WRF and GCM TOE Based On " + gcm_Data_Type]
        avg_WRF_TOE_SDV = avg_Data.variables["SDV of Average TOE of wrf models based on " + wrf_Data_Type][:]
        avg_GCM_TOE_SDV = avg_Data.variables["SDV of Average TOE of gcm models based on " + gcm_Data_Type][:]
        avg_WRF_TOE = avg_WRF_TOE[4:-4,4:-4]
        avg_GCM_TOE = avg_GCM_TOE[4:-4,4:-4]
        avg_WRF_TOE_SDV = avg_WRF_TOE_SDV[4:-4,4:-4]
        avg_GCM_TOE_SDV = avg_GCM_TOE_SDV[4:-4,4:-4]
        avg_Diff = avg_Diff[4:-4,4:-4]
        avg_lon = avg_lon[4:-4,4:-4]
        avg_lat = avg_lat[4:-4,4:-4]
        plt.figure(figsize=(7.5,5))
        
        WRFplot(avg_WRF_TOE, avg_lat, avg_lon,  2010,2050,
                "Average TOE of wrf models based on " + wrf_Data_Type,"Time Of Emergence (years)", "RdYlBu_r")
        plt.savefig(average_Plots_Folder + "//" + "Average TOE of wrf models based on " + wrf_Data_Type + ".png")
        plt.figure(figsize=(7.5,5))
        WRFplot(avg_GCM_TOE, avg_lat, avg_lon,  2010,2050,
                "Average TOE of gcm models based on " + gcm_Data_Type,"Time Of Emergence (years)", "RdYlBu_r")
        plt.savefig(average_Plots_Folder + "//" + "Average TOE of gcm models based on " + gcm_Data_Type + ".png")
        
        plt.figure(figsize=(7.5,5))
        WRFplot(avg_Diff, avg_lat, avg_lon,  -1*amax(avg_Diff), amax(avg_Diff),
                "Difference Between Average WRF and GCM TOE Based On " + gcm_Data_Type, "Time Of Emergence (years)", "RdYlBu_r")
        plt.savefig(average_Plots_Folder + "//" + "Difference Between Average WRF and GCM TOE Based On " + gcm_Data_Type+ ".png")
        
        plt.figure(figsize=(7.5,5))
        
        WRFplot(avg_WRF_TOE_SDV, avg_lat, avg_lon,  -1*amax(avg_WRF_TOE_SDV),amax(avg_WRF_TOE_SDV),
                "SDV of Average TOE of wrf models based on " + wrf_Data_Type,"Standard Deviation", "RdYlBu_r")
        plt.savefig(average_Plots_Folder + "//" + "SDV of Average TOE of wrf models based on " + wrf_Data_Type + ".png")
        
        plt.figure(figsize=(7.5,5))
        
        WRFplot(avg_GCM_TOE_SDV, avg_lat, avg_lon,  -1*amax(avg_GCM_TOE_SDV),amax(avg_GCM_TOE_SDV),
                "SDV of Average TOE of gcm models based on " + gcm_Data_Type,"Standard Deviation", "RdYlBu_r")
        plt.savefig(average_Plots_Folder + "//" + "SDV of Average TOE of gcm models based on " + gcm_Data_Type + ".png")
        avg_Data.close()   
#execute
#define netcdf to store average data
if(os.path.exists("Netcdf_Files" + "\\" + "Average_of_models.nc") == False):
    avg_File = "Netcdf_Files" + "\\" + "Average_of_models.nc"
    ds=clone("Netcdf_Files" + "\\" + "wrf_grid.nc", avg_File)

avg_File = "Netcdf_Files" + "\\" + "Average_of_models.nc"
wrf_Folder = r"Netcdf_Files" + "\\" +"wrf_Netcdf_Files"
gcm_Folder = r"Netcdf_Files" + "\\" +"gcm_Netcdf_Files"




# p
plot_Averages(avg_File)