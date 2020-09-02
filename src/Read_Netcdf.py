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
#from Global_Vars import mean as mn
#from Global_Vars import std 
import matplotlib.pyplot as plt 
import netCDF4 as nc
from netCDF4 import Dataset 
# https://unidata.github.io/netcdf4-python/netCDF4/index.html

# open netcdf file






              

# Get the data from the netcdf file
#
# Note that using these two forms gives different results:
# t90 = ncFile.variables["T90"]
# t90 = ncFile.variables["T90"][:]
#
# The first form creates an object with all the metadata that goes with the variable t90
# The [:] form forces strips the metadata and leaves only the bare numpy array.
# Not sure this is a documented feature, but everyone knows it...
Data_set1 = Dataset("ccsm4_2006-2099_TX90th.nc","r+",format="NETCDF4")
t90 = Data_set1.variables["T90"][:] - 273.15 # convert to deg-C
lat = Data_set1.variables["XLAT"][:]
lon = Data_set1.variables["XLONG"][:]


# get dimensions
ntime = shape(t90)[0] # year index
ny = shape(t90)[1] # south-north - 123
nx = shape(t90)[2] # west-east - 162

# make an array with the year values
year0=2006
year = arange(ntime)+year0


import pylab as lb
if(Data_set1.variables.__str__().find("regressionValues_Slope")==-1):
    regressionValues_Slopes = Data_set1.createVariable("regressionValues_Slope", "float32" , ("south_north","west_east",))
    regressionValues_Slopes.units = "Unkown"
    regressionValues_Yints =  Data_set1.createVariable("regressionValues_Yint", "float32", ("south_north","west_east",))
    regressionValues_Yints.units = "Unknown"


    lb.figure()
#iterate throught the grid points
    Slope=Data_set1.variables["regressionValues_Slope"]
    Yint=Data_set1.variables["regressionValues_Yint"]
    for k in range(ny):
        for i in range(nx):
            #make two arrays to find a linear regression for a single grid point
            x=year
            y=t90[:,[k],[i]]
            Linear_Regression= np.polyfit(x,y,1)
            Linear_Regression = np.squeeze(Linear_Regression)
            x_lin_reg = range(int(year[0]), int(year[ntime-1]))
            predict=np.poly1d(Linear_Regression)
            y_lin_reg = predict(x_lin_reg)
            #save regressions into file 
            Slope[[k],[i]] = Linear_Regression[0]
            Yint[[k],[i]] = Linear_Regression[1]  
            #Plot the linear regression
            lb.plot(x_lin_reg, y_lin_reg, c = 'black')
            #plot the data at that grid point over time
            lb.plot(year,t90[:,[k],[i]])

lb.title("Annual 90th Percentile Tmax")
lb.ylabel("Tmax (deg-C)")
lb.xlabel("Year")
lb.savefig("t90.png")
lb.close()
Data_set1.close()
Data_set2= Dataset("ccsm4_1970-2006_TX90th.nc", "r+", format="NETCDF4")
t90 = Data_set2.variables["T90"][:] - 273.15 # convert to deg-C


import pylab as lb
if(Data_set2.variables.__str__().find("Standard_Deviation_Plus_Mean")==-1):
    stdhigh = Data_set2.createVariable("Standard_Deviation_Plus_Mean", "float32" , ("south_north","west_east",))
    stdhigh.units = "Unkown"
    

    #iterate throught the grid points
    stdhighs=Data_set2.variables["Standard_Deviation_Plus_Mean"]
    for k in range(ny):
        for i in range(nx):
            #find the mean and standard deviations
            Standard_Deviation=np.std(t90[:,[k],[i]])
            mn=np.mean(t90[:,[k],[i]])
            #get mean plus std and record it in file
            stdhighs[[k],[i]]=Standard_Deviation+mn
#close file
Data_set2.close()
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
    


#define ToE NETCDF
fn = "Time_of_Emergence.nc"
#clone historical for the grid data to ToE
ds=clone("ccsm4_1970-2006_TX90th.nc","Time_of_Emergence.nc")

#prepare data retrieval
Data_set4 = Dataset("ccsm4_1970-2006_TX90th.nc","r",format="NETCDF4")
Data_set3 = Dataset("ccsm4_2006-2099_TX90th.nc","r",format="NETCDF4")
#define ToE variable at each grid point
ToE = ds.createVariable("ToE", "float32", ("south_north", "west_east",))
ToE.units = "years"

ToEs=ds.variables["ToE"]


for k in range(ny):
    for i in range(nx):
        slope=Data_set3.variables["regressionValues_Slope"][k][i]
        stdmn=Data_set4.variables["Standard_Deviation_Plus_Mean"][k][i]
        Yint=Data_set3.variables["regressionValues_Yint"][k][i]
        if(slope==0):
            ToE[[k][i]]=3000
        else:
            Time_of_Emergence = (stdmn-Yint)/slope
            ToEs[[k],[i]] = Time_of_Emergence
        
ds.close()
# #to dos
#     #calculate the ToE of Net_cdf at every grid point
#     #visualize ToE at every grid point with basemap library


