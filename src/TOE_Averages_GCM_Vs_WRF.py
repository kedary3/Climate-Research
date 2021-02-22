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


#to compare different wrf and gcm models, we need the toe values on the same grid. 
#To do this we need to create an interpolated toe variable that has the size of the greatest lat and lon combo in the set.
#methodology here is to only upsample the toe data. 

def get_Biggest_GCM_Grid(folder,data_type):
    Max_lat = 0
    Max_lon = 0
    
    biggest_Bounds = 0
    t_Data_Types = ["tasmax90"]
    p_Data_Types = ["pr95"]
    if (data_type in p_Data_Types):
        for file in os.listdir(folder):
            
            gcm_head, gcm_tail = os.path.split(file)
            if(gcm_tail.__str__().find("pr_extr")>=0):
                
                ds = nc.Dataset("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + gcm_tail,"r", format="NETCDF4")
                lats = ds.variables["lat"][:]
                lat = shape(lats)[0]
                lons = ds.variables["lon"][:]
                lon = shape(lons)[0]
                if (lat >= Max_lat and lon >=Max_lon):
                    Max_lat = lat
                    Max_lon = lon
                    biggest_Bounds = gcm_tail 
                
                ds.close()
        return biggest_Bounds
    elif (data_type in t_Data_Types):
        for file in os.listdir(folder):
            gcm_head, gcm_tail = os.path.split(file)
            if(gcm_tail.__str__().find("tasmax_extr")>=0):
                ds = nc.Dataset("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + gcm_tail,"r", format="NETCDF4")
                lats = ds.variables["lat"][:]
                lat = shape(lats)[0]
                lons = ds.variables["lon"][:]
                lon = shape(lons)[0]
                if (lat >= Max_lat and lon >=Max_lon):
                    biggest_Bounds = gcm_tail 
                
                ds.close()
        return biggest_Bounds
    return "invalid data_type"

# def get_Biggest_WRF_Grid(folder,data_type):
#     Max_lat = 0
#     Max_lon = 0
    
#     biggest_Bounds = 0
#     t_Data_Types = ["T2MAX90"]
#     p_Data_Types = ["PREC95"]
#     if (data_type in p_Data_Types):
#         for file in os.listdir(folder):
#             wrf_head, wrf_tail = os.path.split(file)
#             if(wrf_tail.__str__().find("PREC_extr")>=0):
#                 ds = nc.Dataset("Netcdf_Files" + "\\" + "wrf_Netcdf_Files" + "\\" + wrf_tail,"r", format="NETCDF4")
#                 lats = ds.variables["XLAT"][:]
#                 lat = shape(lats)[0]
#                 lons = ds.variables["XLONG"][:]
#                 lon = shape(lons)[0]
#                 if (lat >= Max_lat and lon >=Max_lon):
#                     biggest_Bounds = wrf_tail 
                
#                 ds.close()
#         return biggest_Bounds
#     elif (data_type in t_Data_Types):
#         for file in os.listdir(folder):
#             wrf_head, wrf_tail = os.path.split(file)
#             if(wrf_tail.__str__().find("T2MAX_extr")>=0):
#                 ds = nc.Dataset("Netcdf_Files" + "\\" + "wrf_Netcdf_Files" + "\\" + wrf_tail,"r", format="NETCDF4")
#                 lats = ds.variables["XLAT"][:]
#                 lat = shape(lats)[0]
#                 lons = ds.variables["XLONG"][:]
#                 lon = shape(lons)[0]
#                 if (lat >= Max_lat and lon >=Max_lon):
#                     biggest_Bounds = wrf_tail 
                
#                 ds.close()
#         return biggest_Bounds
#     return "invalid data_type"
def interpolate_GCM_ToE(folder,data_type):
    
    gcm_t_Data_Types = ["tasmax90", "tasmaxx"]
    gcm_p_Data_Types = ["pr95", "prx"]
    for file in os.listdir(folder):
        gcm_head, gcm_tail = os.path.split(file)
        if(gcm_tail.__str__().find("pr_extr")>=0 and data_type in gcm_p_Data_Types):
                
            Idata=Dataset("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + get_Biggest_GCM_Grid(folder,data_type), "r", format="NETCDF4")
            
            inter_lon=Idata.variables["lon"][:] 
            inter_lat=Idata.variables["lat"][:] 
            data=Dataset("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + gcm_tail,"r+", format="NETCDF4")
            if(data.dimensions.__str__().find("big-lat2")==-1):
                data.createDimension( "big-lon2" , size=inter_lon.shape[0])
                data.createDimension( "big-lat2" , size=inter_lat.shape[0])
            if(data.variables.__str__().find("GCM_Interpolation_of_ToE_based_on_"+data_type)<0):
                inter_ToE = data.createVariable("GCM_Interpolation_of_ToE_based_on_"+data_type,
                                          "float32", ("big-lat2", "big-lon2",))
                inter_ToE.units = "years"
            ToE=data.variables["ToE_for_"+data_type][:]
            lon=data.variables["lon"][:] 
            lat=data.variables["lat"][:] 
            interpolated_TOE = data.variables["GCM_Interpolation_of_ToE_based_on_"+data_type][:]
            
            # Here use x as lon and y as lat and Z as data field from GCM
            x = lon
            y = lat
            z = ToE[:,:]
            
            # here X2 and Y2 (upper case) are lon and lat from WRF grid.nc
            # note that lat and lon are 2d arrays 
            #x2 = linspace(0, 1, 123)
            #y2 = linspace(0, 1, 162)
            #X2,Y2 = meshgrid(x2,y2) # from grid.nc
            
            
            
            interp = interp2d(lon, lat, z, kind='cubic')
            # Since X2 and Y2 are 2-d arrays for irregular grid need to do point by point
            #    **THERE MAY BE A BETTER WAY?
            inter_ToE=zeros_like(1,shape=(inter_lon.shape[0],inter_lat.shape[0]))
            
            for i in range(inter_lon.shape[0]):
                for j in range(inter_lat.shape[0]):
                    inter_ToE[i,j]=interp(inter_lon[i],inter_lat[j])
                    print(inter_ToE[i,j])
                    interpolated_TOE[[j],[i]] = inter_ToE[i,j]
            data.close()
            Idata.close()
        elif(gcm_tail.__str__().find("tasmax_extr")>=0 and data_type in gcm_t_Data_Types):
            Idata=Dataset("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + get_Biggest_GCM_Grid(folder,data_type), "r", format="NETCDF4")
            
            inter_lon=Idata.variables["lon"][:] 
            inter_lat=Idata.variables["lat"][:] 
            data=Dataset("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + gcm_tail,"r+", format="NETCDF4")
            if(data.dimensions.__str__().find("big-lat2")==-1):
                data.createDimension( "big-lon2" , size=inter_lon.shape[0])
                data.createDimension( "big-lat2" , size=inter_lat.shape[0])
            if(data.variables.__str__().find("GCM_Interpolation_of_ToE_based_on_"+data_type)<0):
                inter_ToE = data.createVariable("GCM_Interpolation_of_ToE_based_on_"+data_type,
                                          "float32", ("big-lat2", "big-lon2",))
                inter_ToE.units = "years"
            ToE=data.variables["ToE_for_"+data_type][:]
            lon=data.variables["lon"][:] 
            lat=data.variables["lat"][:] 
            interpolated_TOE = data.variables["GCM_Interpolation_of_ToE_based_on_"+data_type][:]
            
            # Here use x as lon and y as lat and Z as data field from GCM
            x = lon
            y = lat
            z = ToE[:,:]
            
            # here X2 and Y2 (upper case) are lon and lat from WRF grid.nc
            # note that lat and lon are 2d arrays 
            #x2 = linspace(0, 1, 123)
            #y2 = linspace(0, 1, 162)
            #X2,Y2 = meshgrid(x2,y2) # from grid.nc
            
            
            
            interp = interp2d(lon, lat, z, kind='cubic')
            # Since X2 and Y2 are 2-d arrays for irregular grid need to do point by point
            #    **THERE MAY BE A BETTER WAY?
            inter_ToE=zeros_like(1,shape=(inter_lon.shape[0],inter_lat.shape[0]))
            
            for i in range(inter_lon.shape[0]):
                for j in range(inter_lat.shape[0]):
                    inter_ToE[i,j]=interp(inter_lon[i],inter_lat[j])
                    interpolated_TOE[[j],[i]] = inter_ToE[i,j]
            data.close()
            Idata.close()
def interpolate_GCM_SDV(folder,data_type):
    
    for file in os.listdir(folder):
        gcm_head, gcm_tail = os.path.split(file)
        data=Dataset("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + gcm_tail,"r+", format="NETCDF4")
        if(data.variables.__str__().find("Interpolation_of_SDV_based_on_"+data_type)<0):
            inter_SDV = data.createVariable("Interpolation_of_SDV_based_on_"+data_type,
                                      "float32", ("south_north", "west_east",))
            
        SDV=data.variables["Standard_Deviations_for_"+data_type][:]
        lon=data.variables["lon"][:] 
        lat=data.variables["lat"][:] 
        interpolated_SDV = data.variables["Interpolation_of_SDV_based_on_"+data_type][:]
        
        # Here use x as lon and y as lat and Z as data field from GCM
        x = lon
        y = lat
        Z = SDV[:,:]
        
        # here X2 and Y2 (upper case) are lon and lat from WRF grid.nc
        # note that lat and lon are 2d arrays 
        #x2 = linspace(0, 1, 123)
        #y2 = linspace(0, 1, 162)
        #X2,Y2 = meshgrid(x2,y2) # from grid.nc
        
        Idata=Dataset("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + get_Biggest_GCM_Grid(folder,data_type), "r", format="NETCDF4")
        
        
        inter_lon=Idata.variables["lon"][:] 
        inter_lat=Idata.variables["lat"][:] 
        
        
        interp = interp2d(lon, lat, z, kind='cubic')
        
        # Since X2 and Y2 are 2-d arrays for irregular grid need to do point by point
        #    **THERE MAY BE A BETTER WAY?
        inter_SDV=zeros_like(inter_lon)
        for i in arange(inter_lon.shape[0]):
            for j in arange(inter_lat.shape[0]):
                inter_SDV[i,j]=interp(inter_lon[i,j],inter_lat[i,j])   
                
                interpolated_SDV[[i],[j]] = inter_SDV[i,j]
        data.close()
        Idata.close()
                               
#method to find the average of all gcm and wrf models' TOE for the PNW.  


def Calculate_Average_GCM_TOE(file, data_type, folder):
    ncFile = nc.Dataset(file,"r+", format="NETCDF4")
    AVG_TOE = ncFile.variables["GCM_Average_TOE_for_"+data_type][:]
    AVG_SDV = ncFile.variables["GCM_Average_SDV_for_"+data_type][:]
    lat = ncFile.variables["lat"][:]
    lon = ncFile.variables["lon"][:]
    interpolate_GCM_ToE(folder,data_type)
    # interpolate_GCM_SDV(folder,data_type)
    
    n=0
    ny = shape(lat)[0] # south-north - 123 - wrf / south-north - 21 - gcm
    nx = shape(lon)[0] # west-east - 162 - wrf / west-east - 33 - gcm
    print(ny,nx)
    SUM_TOE = np.zeros((ny,nx))
    SUM_SDV = np.zeros((ny,nx))
    for netcdf_file in os.listdir(folder):
        head, tail = os.path.split(netcdf_file)
        if (tail.find("pr_extr")>=0):
            n+=1
            print(tail)
            ds = nc.Dataset("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + tail, "r", format = "NETCDF4")
            toe = ds.variables["GCM_Interpolation_of_ToE_based_on_"+data_type][:]
            # sdv = ds.variables["Standard_Deviations_for_"+data_type][:]
            for k in range(ny):
                for i in range(nx):
                    SUM_TOE[k][i] += toe[k][i] 
                    # SUM_SDV[k][i] += sdv[k][i] 
            ds.close() 
        
    if n == 0: 
        print("div by zero")
        return 1
    for k in range(ny):
        for i in range(nx):
            AVG_TOE[k][i] = SUM_TOE[k][i]/n
            # AVG_SDV[k][i] = SUM_SDV[k][i]/n
    ncFile.close()
    
def Calculate_Average_WRF_TOE(file, data_type, folder):
    ncFile = nc.Dataset(file,"r+", format="NETCDF4")
    AVG_TOE = ncFile.variables["WRF_Average_TOE_for_"+data_type]
    AVG_SDV = ncFile.variables["WRF_Average_SDV_for_"+data_type]
    lats = ncFile.variables["XLAT"][:]
    lons = ncFile.variables["XLONG"][:]
    
    n=np.size(os.listdir(folder))
    if n == 0: 
        print("div by zero")
        return 1
    ny = shape(lats)[0] # south-north - 123 - wrf / south-north - 21 - gcm
    nx = shape(lons)[1] # west-east - 162 - wrf / west-east - 33 - gcm
    print(ny,nx)
    SUM_TOE = np.zeros_like((nx,ny))
    SUM_SDV = np.zeros_like((nx,ny))
    for netcdf_file in os.listdir(folderr):
        head, tail = os.path.split(netcdf_file)
        ds = nc.Dataset(netcdf_file, "r", format = "NETCDF4")
        toe = ds.variables["ToE_for_"+data_type][:]
        sdv = ds.variables["Standard_Deviations_for_",data_type][:]
        for k in range(ny):
            for i in range(nx):
                SUM_TOE[k][i] += toe[k][i] 
                SUM_SDV[k][i] += sdv[k][i] 
        ds.close()  
    for k in range(ny):
        for i in range(nx):
            AVG_TOE[[k],[i]] = SUM_TOE[k][i]/n
            AVG_SDV[[k],[i]] = SUM_SDV[k][i]/n
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

def GCM_clone(src_file, trg_file): #function to copy attributes, variables, and dimensions from one NETCDF to another
    src = nc.Dataset(src_file)
    trg = nc.Dataset(trg_file, mode='w')

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        # if(var.__str__().find("lat")>=0 or var.__str__().find("lon")>=0 or 
        #    var.__str__().find("lat_bnds")>=0 or var.__str__().find("lon_bnds")>=0 or
        #    var.__str__().find("time")>=0 or var.__str__().find("time_bnds")>=0):
            trg.createVariable(name, var.dtype, var.dimensions)

            # Copy the variable attributes
            trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})
    
            # Copy the variables values (as 'f4' eventually)
            trg.variables[name][:] = src.variables[name][:]
    # return the file and close the source
    src.close
    return trg

def main():
    #define list of gcm,wrf files
    wrf_Folder = r"Netcdf_Files" + "\\" +"wrf_Netcdf_Files"
    gcm_Folder = r"Netcdf_Files" + "\\" +"gcm_Netcdf_Files"
    
    
    #define list of data types
    GCM_t_Data_Types = ["tasmax90"]
    GCM_p_Data_Types = ["pr95"] 
    WRF_t_Data_Types = ["T2MAX90"]
    WRF_p_Data_Types = ["PREC95"] 
    
    
    #for each wrf and gcm file lists, generate toe Average data for each data type
    #define netcdf files to store average data
    #GCM file creation
    if(os.path.exists("Netcdf_Files" + "\\" + "Average_TOE_GCM.nc") == False):
        AV_GCM_File = "Netcdf_Files" + "\\" + "Average_TOE_GCM.nc"
        fs=GCM_clone("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + get_Biggest_GCM_Grid(gcm_Folder,"tasmax90"), AV_GCM_File)

    AV_GCM_File = "Netcdf_Files" + "\\" + "Average_TOE_GCM.nc"
    ds = nc.Dataset(AV_GCM_File, "r+", format = "NETCDF4")
    for GCM_p_Data_Type in GCM_p_Data_Types:
        if(ds.variables.__str__().find("GCM_Average_TOE_for_"+GCM_p_Data_Type)==-1):
            AVToE = ds.createVariable("GCM_Average_TOE_for_"+GCM_p_Data_Type, "float32", ("lat","lon",))
            AVSDV = ds.createVariable("GCM_Average_SDV_for_"+GCM_p_Data_Type, "float32", ("lat","lon",))
        
        ds.close()
        
        Calculate_Average_GCM_TOE("Netcdf_Files" + "\\" + "Average_TOE_GCM.nc", GCM_p_Data_Type, gcm_Folder)
        
                    
        
    # for GCM_t_Data_Type in GCM_t_Data_Types:
    #     if(AV_GCM_File.variables.__str__().find("GCM_Average_TOE_for"+GCM_t_Data_Type)==-1):
    #         ToE = AV_GCM_File.createVariable("GCM_Average_TOE_for"+GCM_t_Data_Type, "float32", ("lat","lon",))
    #         ToE.units = "years"
    #     else:
    #             print("Average GCM TOE for",GCM_t_Data_Type,  "data already generated")
    #     Calculate_Average_TOE("Netcdf_Files" + "\\" + "Average_TOE_GCM.nc",
    #                           GCM_t_Data_Type, gcm_Folder)
    # #WRF
    
    # if(os.path.exists("Netcdf_Files" + "\\" + "Average_TOE_WRF.nc") == False):  
    #     AV_WRF_File = "Netcdf_Files" + "\\" + "Average_TOE_WRF.nc"
    #     ds=WRF_clone(os.listdir(wrf_Folder)[0], AV_WRF_File)
    # AV_WRF_File = Dataset(ds, "r+", format = "NETCDF4")
    # for WRF_p_Data_Type in WRF_p_Data_Types:
    #     if(AV_WRF_File.variables.__str__().find("WRF_Average_TOE_for"+WRF_p_Data_Type)==-1):
    #         ToE = AV_WRF_File.createVariable("WRF_Average_TOE_for"+WRF_p_Data_Type, "float32", ("XLAT","XLONG",))
    #         ToE.units = "years"
    #     else:
    #             print("Average WRF TOE for",WRF_p_Data_Type,  "data already generated")
    #     Calculate_Average_TOE("Netcdf_Files" + "\\" + "Average_TOE_WRF.nc",
    #                           WRF_p_Data_Type, wrf_Folder)
    # for WRF_t_Data_Type in WRF_t_Data_Types:
    #     if(AV_WRF_File.variables.__str__().find("WRF_Average_TOE_for"+WRF_t_Data_Type)==-1):
    #         ToE = AV_WRF_File.createVariable("WRF_Average_TOE_for"+WRF_t_Data_Type, "float32", ("XLAT","XLONG",))
    #         ToE.units = "years"
    #     else:
    #             print("Average WRF TOE for",WRF_t_Data_Type,  "data already generated")
    #     Calculate_Average_TOE("Netcdf_Files" + "\\" + "Average_TOE_WRF.nc",
    #                           WRF_t_Data_Type, wrf_Folder)
#execute
main()