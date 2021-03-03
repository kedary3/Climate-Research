# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 06:00:42 2021

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




#define list of data types
GCM_t_Data_Types = ["tasmax90"]
GCM_p_Data_Types = ["pr95"] 
WRF_t_Data_Types = ["T2MAX90"]
WRF_p_Data_Types = ["PREC95"] 
      
#method to find the average of all gcm and wrf models' TOE for the PNW.  


def Calculate_Average_GCM_TOE(file, data_type, folder):
    ncFile = nc.Dataset(file,"r+", format="NETCDF4")
    AVG_TOE = ncFile.variables["GCM_Average_TOE_for_"+data_type]
    AVG_SDV = ncFile.variables["GCM_Average_SDV_for_"+data_type]
    
    
    
    n=0
    ny = 123 # south-north - 123 - wrf / south-north - 21 - gcm
    nx = 162 # west-east - 162 - wrf / west-east - 33 - gcm
    SUM_TOE = np.zeros_like(2,shape=(ny,nx), dtype= float)
    SUM_SDV = np.zeros_like(2,shape=(ny,nx), dtype= float)
    for netcdf_file in os.listdir(folder):
        head, tail = os.path.split(netcdf_file)
        if (data_type in GCM_p_Data_Types and tail.find("pr_extr")>=0):
            n+=1
            ds = nc.Dataset("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + tail, "r", format = "NETCDF4")
            toe = ds.variables["Interpolated ToE data based on " + data_type][:]
            sdv = ds.variables["Interpolated SDV data based on " + data_type][:]
            for k in range(ny):
                for i in range(nx):
                    
                    SUM_TOE[[k],[i]] += toe[[k],[i]] 
                    SUM_SDV[[k],[i]] += sdv[[k],[i]] 
            ds.close() 
        if (data_type in GCM_t_Data_Types and tail.find("tasmax_extr")>=0):
            n+=1
            ds = nc.Dataset("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + tail, "r", format = "NETCDF4")
            toe = ds.variables["Interpolated ToE data based on " + data_type][:]
            sdv = ds.variables["Interpolated SDV data based on " + data_type][:]
            for k in range(ny):
                for i in range(nx):
                    
                    SUM_TOE[[k],[i]] += toe[[k],[i]] 
                    SUM_SDV[[k],[i]] += sdv[[k],[i]] 
            ds.close() 
    if n == 0: 
        print("div by zero")
        return 1
    for k in range(ny):
        for i in range(nx):
            AVG_TOE[[k],[i]] = np.float32(SUM_TOE[[k],[i]][0]/n)
            AVG_SDV[[k],[i]] = np.float32(SUM_SDV[[k],[i]][0]/n)
    ncFile.close()
    
def Calculate_Average_WRF_TOE(file, data_type, folder):
    ncFile = nc.Dataset(file,"r+", format="NETCDF4")
    AVG_TOE = ncFile.variables["WRF_Average_TOE_for_"+data_type]
    AVG_SDV = ncFile.variables["WRF_Average_SDV_for_"+data_type]
    lats = ncFile.variables["XLAT"][:]
    lons = ncFile.variables["XLONG"][:]
    
    
    ny = shape(lats)[0] # south-north - 123 - wrf / south-north - 21 - gcm
    nx = shape(lons)[1] # west-east - 162 - wrf / west-east - 33 - gcm
    
    SUM_TOE = np.zeros_like(2,shape=(ny,nx), dtype= float)
    SUM_SDV = np.zeros_like(2,shape=(ny,nx), dtype= float)
    n = 0
    for netcdf_file in os.listdir(folder):
        head, tail = os.path.split(netcdf_file)
        if (data_type in WRF_p_Data_Types and tail.find("PREC_extr")>=0):
            n+=1
            ds = nc.Dataset("Netcdf_Files" + "\\" + "wrf_Netcdf_Files" + "\\" + tail, "r", format = "NETCDF4")
            toe = ds.variables["ToE_for_"+data_type][:]
            sdv = ds.variables["Standard_Deviations_for_"+data_type][:]
            for k in range(ny):
                for i in range(nx):
                    SUM_TOE[[k],[i]] += toe[[k],[i]] 
                    SUM_SDV[[k],[i]] += sdv[[k],[i]] 
            ds.close()
        if (data_type in WRF_t_Data_Types and tail.find("T2MAX_extr")>=0):
            n+=1
            ds = nc.Dataset("Netcdf_Files" + "\\" + "wrf_Netcdf_Files" + "\\" + tail, "r", format = "NETCDF4")
            toe = ds.variables["ToE_for_"+data_type][:]
            sdv = ds.variables["Standard_Deviations_for_"+data_type][:]
            for k in range(ny):
                for i in range(nx):
                    SUM_TOE[[k],[i]] += toe[[k],[i]] 
                    SUM_SDV[[k],[i]] += sdv[[k],[i]] 
            ds.close()  
    for k in range(ny):
        for i in range(nx):
            AVG_TOE[[k],[i]] = SUM_TOE[[k],[i]]/n
            AVG_SDV[[k],[i]] = SUM_SDV[[k],[i]]/n
    ncFile.close()
# define clone function to match netcdf information
def WRF_clone(src_file, trg_file): #function to copy attributes, variables, and dimensions from one NETCDF to another
    src = nc.Dataset(src_file)
    trg = nc.Dataset(trg_file, mode='w')

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        # if(var.__str__().find("XLAT")>=0 or var.__str__().find("XLONG")>=0 or 
        #    var.__str__().find("XTIME")>=0 or var.__str__().find("XTIME_bnds")>=0):
            trg.createVariable(name, var.dtype, var.dimensions)

            # Copy the variable attributes
            trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})
    
            # Copy the variables values (as 'f4' eventually)
            trg.variables[name][:] = src.variables[name][:]
    # return the file and close the source
    src.close
    return trg


def main(grid_file):
    #define list of gcm,wrf files
    wrf_Folder = r"Netcdf_Files" + "\\" +"wrf_Netcdf_Files"
    gcm_Folder = r"Netcdf_Files" + "\\" +"gcm_Netcdf_Files"
    
    
    
    
    
    #for each wrf and gcm file lists, generate toe Average data for each data type
    #define netcdf file to store average data
    #GCM file creation
    if(os.path.exists("Netcdf_Files" + "\\" + "Average_TOE.nc") == False):
        AV_TOE_File = "Netcdf_Files" + "\\" + "Average_TOE.nc"
        fs=WRF_clone(grid_file, AV_TOE_File)

    AV_File = "Netcdf_Files" + "\\" + "Average_TOE.nc"
    ds = nc.Dataset(AV_File, "r+", format = "NETCDF4")
    for GCM_p_Data_Type in GCM_p_Data_Types:
        if(ds.variables.__str__().find("GCM_Average_TOE_for_"+GCM_p_Data_Type)==-1):
            AVToE = ds.createVariable("GCM_Average_TOE_for_"+GCM_p_Data_Type, "float32", ("south_north","west_east",))
            AVSDV = ds.createVariable("GCM_Average_SDV_for_"+GCM_p_Data_Type, "float32", ("south_north","west_east",))
        else:
            print("Average GCM TOE for",GCM_p_Data_Type,  "data already generated")
        ds.close()
        Calculate_Average_GCM_TOE("Netcdf_Files" + "\\" + "Average_TOE.nc",
                                  GCM_p_Data_Type, gcm_Folder)
    
                    
    ds = nc.Dataset(AV_File, "r+", format = "NETCDF4")   
    for GCM_t_Data_Type in GCM_t_Data_Types:
        ds = nc.Dataset(AV_File, "r+", format = "NETCDF4")
        if(ds.variables.__str__().find("GCM_Average_TOE_for_"+GCM_t_Data_Type)==-1):
            AVToE = ds.createVariable("GCM_Average_TOE_for_"+GCM_t_Data_Type, "float32", ("south_north","west_east",))
            AVSDV = ds.createVariable("GCM_Average_SDV_for_"+GCM_t_Data_Type, "float32", ("south_north","west_east",))
        else:
            print("Average GCM TOE for",GCM_t_Data_Type,  "data already generated")
        ds.close()
        Calculate_Average_GCM_TOE("Netcdf_Files" + "\\" + "Average_TOE.nc",
                              GCM_t_Data_Type, gcm_Folder)
    
    #WRF
    
    ds = nc.Dataset(AV_File, "r+", format = "NETCDF4")
    for WRF_p_Data_Type in WRF_p_Data_Types:
        if(ds.variables.__str__().find("WRF_Average_TOE_for_"+WRF_p_Data_Type)==-1):
            ToE = ds.createVariable("WRF_Average_TOE_for_"+WRF_p_Data_Type, "float32", ("south_north","west_east",))
            ToE.units = "years"
            AVSDV = ds.createVariable("WRF_Average_SDV_for_"+WRF_p_Data_Type, "float32", ("south_north","west_east",))
        else:
            print("Average WRF TOE for",WRF_p_Data_Type,  "data already generated")
        ds.close()
        Calculate_Average_WRF_TOE("Netcdf_Files" + "\\" + "Average_TOE.nc",
                              WRF_p_Data_Type, wrf_Folder)
    ds = nc.Dataset(AV_File, "r+", format = "NETCDF4")
    for WRF_t_Data_Type in WRF_t_Data_Types:
        if(ds.variables.__str__().find("WRF_Average_TOE_for_"+WRF_t_Data_Type)==-1):
            ToE = ds.createVariable("WRF_Average_TOE_for_"+WRF_t_Data_Type, "float32", ("south_north","west_east",))
            ToE.units = "years"
            AVSDV = ds.createVariable("WRF_Average_SDV_for_"+WRF_t_Data_Type, "float32", ("south_north","west_east",))
        else:
            print("Average WRF TOE for",WRF_t_Data_Type,  "data already generated")
        ds.close()
        Calculate_Average_WRF_TOE("Netcdf_Files" + "\\" + "Average_TOE.nc",
                              WRF_t_Data_Type, wrf_Folder)
    
#execute
grid_file = "Netcdf_files\\wrf_grid.nc"
main(grid_file)