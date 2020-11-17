# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 11:22:28 2020

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
import pylab as lb
from numpy import *
from scipy.interpolate import interp2d

def generate_GCM_PREC_ToE_Data(file, Precipitation_Type):
    
    #check input before working on files
    
    if(Precipitation_Type !="prx" and Precipitation_Type !="pr95"):
        print("invalid Precipitation data type")
        return 1 
    
    #open netcdf file
    
    ncFile = Dataset(file,"r+", format="NETCDF4")
    
        
    # get the data to plot
    PREC=ncFile.variables[Precipitation_Type][:]*24*60*60 # mm/day
    lats = ncFile.variables["lat"][:]
    lons = ncFile.variables["lon"][:]
    
    #get time axis and grid values
    year=shape(PREC)[0]
    year0=1970
    year=arange(year)+year0
    year1=1999
    year2=2099
    ny = shape(PREC)[1] # south-north - 123
    nx = shape(PREC)[2] # west-east - 162
    
    #Prepare Variables for ToE Calculation
    # if(ncFile.variables.__str__().find("regressionValues_Slope_for_" + Precipitation_Type)==-1):
    #     regressionValues_Slopes = ncFile.createVariable("regressionValues_Slope_for_" + Precipitation_Type,
    #                                                     "float32" , ("lat","lon",))
    #     regressionValues_Slopes.units = "kg per m^2 per year"
    #     regressionValues_Yints =  ncFile.createVariable("regressionValues_Yint_for_" + Precipitation_Type,
    #                                                     "float32", ("lat","lon",))
    #     regressionValues_Yints.units = "kg per m^2"
    # if(ncFile.variables.__str__().find("Standard_Deviations_for_" + Precipitation_Type)==-1):
    #     std = ncFile.createVariable("Standard_Deviations_for_" + Precipitation_Type, "float32" , ("lat","lon",))
    #     std.units = "kg per m^2"
    # if(ncFile.variables.__str__().find("Means_for_" + Precipitation_Type)==-1):
    #     mean = ncFile.createVariable("Means_for_" + Precipitation_Type, "float32" , ("lat","lon",))
    #     mean.units = "kg per m^2" 
    # #define ToE variable at each grid point
    # if(ncFile.variables.__str__().find("ToE_for_" + Precipitation_Type)==-1):
    #     ToE = ncFile.createVariable("ToE_for_" + Precipitation_Type, "float32", ("lat","lon",))
    #     ToE.units = "years"
    # else:
    #     print("ToE data already generated for " + Precipitation_Type)
    #     return 1
    #iterate throught the grid points and assign the new variables their respective values
    
    Slopes   = ncFile.variables["regressionValues_Slope_for_" + Precipitation_Type]
    Yints    = ncFile.variables["regressionValues_Yint_for_" + Precipitation_Type]
    stds     = ncFile.variables["Standard_Deviations_for_" + Precipitation_Type]
    means    = ncFile.variables["Means_for_" + Precipitation_Type]
    
    lb.figure()
    for k in range(ny):
        for i in range(nx):
            
            #find the mean and standard deviations
            Standard_Deviation=np.std(PREC[:,[k],[i]][0:31])
            mn=np.mean(PREC[:,[k],[i]][0:31])
            
            #get mean plus std and record it in file
            means[[k],[i]]=mn
            stds[[k],[i]]=Standard_Deviation
            
            #make two arrays to find a linear regression for a single grid point
            x = year[30:130]
            y = PREC[:,[k],[i]][30:130]
            Linear_Regression = np.polyfit(x,y,1)
            Linear_Regression = np.squeeze(Linear_Regression)
            x_lin_reg = range(int(year[30]), int(year[129]))
            predict=np.poly1d(Linear_Regression)
            y_lin_reg = predict(x_lin_reg)
            #save regressions into file 
            Slopes[[k],[i]] = Linear_Regression[0]
            Yints[[k],[i]] = Linear_Regression[1]  
            slope = Linear_Regression[0]
            Yint  = Linear_Regression[1]  
            #Plot the linear regression
            lb.plot(x_lin_reg, y_lin_reg, c = 'black')
            #plot the data at that grid point over time
            lb.plot(year,PREC[:,[k],[i]])
            
            
    lb.title("prx from " + file.__str__())
    lb.ylabel("kg/s*m^2")
    lb.xlabel("Year")
    lb.show()       
    #close file
    ncFile.close()
    
import os
folder = r'Netcdf_Files'
#for each wrf and gcm file, generate toe data for each data type and get interpolated gcm data
for file in os.scandir(folder):
    head, tail = os.path.split(file) #tail gives the file name and type
    if(file.__str__().find("prextr")!=-1): #if file has precipitation data 
        generate_GCM_PREC_ToE_Data("Netcdf_Files" + "\\" + tail, "prx")