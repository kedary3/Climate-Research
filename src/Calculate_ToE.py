# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 00:44:39 2020

@author: Kedar
"""


#todo
    #generate plots of  at a specific gridpoints
        
    #change year limit around
    #median ToEs
    #change PRECIP calculation to appendix 4 calculation
        #compare median gcm vs median wrf
    #take out border of toe plots for masking precipitation
    #take out border of delta plots for masking  

#ASKS
    #ask about ddof in std

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


def generate_GCM_TASMAX_ToE_Data(file, Temperature_Type):
    
    #check input before working on files
    
    if(Temperature_Type !="tasmaxx" and Temperature_Type !="tasmaxn" and Temperature_Type != "tasmax90"):
        print("invalid Temperature data type")
        return 1 
    
    #open netcdf file
    
    ncFile = Dataset(file,"r+", format="NETCDF4")
    
        
    # get the data to plot
    TEMP=ncFile.variables[Temperature_Type][:]-273.15 #convert to Kelvin
    lats = ncFile.variables["lat"][:]
    lons = ncFile.variables["lon"][:]
    
    #get time axis and grid values
    year=shape(TEMP)[0]
    year0=1970
    year=arange(year)+year0
    year1=1999
    year2=2099
    yearLimit=2300
    ny = shape(TEMP)[1] # south-north - 123 - wrf / south-north - 21 - gcm
    nx = shape(TEMP)[2] # west-east - 162 - wrf / west-east - 33 - gcm
    #splice the variable array along the time axis, separating model and historical data
    #31 is chosen because python is exclusive on the second range number
    Historical_TEMP = TEMP[:,1,1][0:31]
    plt.plot(year[0:31],Historical_TEMP)
    Model_TEMP = TEMP[:,1,1][30:130]
    plt.plot(year[30:130],Model_TEMP)
    plt.show()
    
    #Prepare Variables for ToE Calculation
    if(ncFile.variables.__str__().find("regressionValues_Slope_for_" + Temperature_Type)==-1):
        regressionValues_Slopes = ncFile.createVariable("regressionValues_Slope_for_" + Temperature_Type,
                                                        "float32" , ("lat","lon",))
        regressionValues_Slopes.units = "kg per m^2 per year"
        regressionValues_Yints =  ncFile.createVariable("regressionValues_Yint_for_" + Temperature_Type,
                                                        "float32", ("lat","lon",))
        regressionValues_Yints.units = "kg per m^2"
    if(ncFile.variables.__str__().find("Standard_Deviations_for_" + Temperature_Type)==-1):
        std = ncFile.createVariable("Standard_Deviations_for_" + Temperature_Type, "float32" , ("lat","lon",))
        std.units = "kg per m^2"
    if(ncFile.variables.__str__().find("Means_for_" + Temperature_Type)==-1):
        mean = ncFile.createVariable("Means_for_" + Temperature_Type, "float32" , ("lat","lon",))
        mean.units = "kg per m^2" 
    #define ToE variable at each grid point
    if(ncFile.variables.__str__().find("ToE_for_" + Temperature_Type)==-1):
        ToE = ncFile.createVariable("ToE_for_" + Temperature_Type, "float32", ("lat","lon",))
        ToE.units = "years"
    else:
        print("ToE data already generated for " + Temperature_Type)
        return 1
    #iterate throught the grid points and assign the new variables their respective values
    
    Slopes   = ncFile.variables["regressionValues_Slope_for_" + Temperature_Type]
    Yints    = ncFile.variables["regressionValues_Yint_for_" + Temperature_Type]
    stds     = ncFile.variables["Standard_Deviations_for_" + Temperature_Type]
    means    = ncFile.variables["Means_for_" + Temperature_Type]
    ToEs     = ncFile.variables["ToE_for_" + Temperature_Type]
    lb.figure()
    for k in range(ny):
        for i in range(nx):
            
            #find the mean and standard deviations
            Standard_Deviation=np.std(TEMP[:,[k],[i]][0:31], ddof=1)
            mn=np.mean(TEMP[:,[k],[i]][0:31])
            
            #get mean plus std and record it in file
            means[[k],[i]]=mn
            stds[[k],[i]]=Standard_Deviation
            
            #make two arrays to find a linear regression for a single grid point
            x = year[30:130]
            y = TEMP[:,[k],[i]][30:130]
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
            
            #check for divByZero
            if(slope == 0):
                ToEs[[k],[i]]=year2
            else:
                #calculate ToE
                Time_of_Emergence = year1 + Standard_Deviation/np.abs(slope)
                if(Time_of_Emergence>year2):
                    ToEs[[k],[i]]=year2
                else:
                    ToEs[[k],[i]] = Time_of_Emergence
            
    lb.title(Temperature_Type)
    lb.ylabel("Tmax (deg-C)")
    lb.xlabel("Year")
    lb.show()       
    #close file
    ncFile.close()

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
    if(ncFile.variables.__str__().find("regressionValues_Slope_for_" + Precipitation_Type)==-1):
        regressionValues_Slopes = ncFile.createVariable("regressionValues_Slope_for_" + Precipitation_Type,
                                                        "float32" , ("lat","lon",))
        regressionValues_Slopes.units = "mm per day per year"
        regressionValues_Yints =  ncFile.createVariable("regressionValues_Yint_for_" + Precipitation_Type,
                                                        "float32", ("lat","lon",))
        regressionValues_Yints.units = "mm/day"
    if(ncFile.variables.__str__().find("Standard_Deviations_for_" + Precipitation_Type)==-1):
        std = ncFile.createVariable("Standard_Deviations_for_" + Precipitation_Type, "float32" , ("lat","lon",))
        std.units = "mm/day"
    if(ncFile.variables.__str__().find("Means_for_" + Precipitation_Type)==-1):
        mean = ncFile.createVariable("Means_for_" + Precipitation_Type, "float32" , ("lat","lon",))
        mean.units = "mm/day" 
    #define ToE variable at each grid point
    if(ncFile.variables.__str__().find("ToE_for_" + Precipitation_Type)==-1):
        ToE = ncFile.createVariable("ToE_for_" + Precipitation_Type, "float32", ("lat","lon",))
        ToE.units = "years"
    else:
        print("ToE data already generated for " + Precipitation_Type)
        return 1
    #iterate throught the grid points and assign the new variables their respective values
    
    Slopes   = ncFile.variables["regressionValues_Slope_for_" + Precipitation_Type]
    Yints    = ncFile.variables["regressionValues_Yint_for_" + Precipitation_Type]
    stds     = ncFile.variables["Standard_Deviations_for_" + Precipitation_Type]
    means    = ncFile.variables["Means_for_" + Precipitation_Type]
    ToEs     = ncFile.variables["ToE_for_" + Precipitation_Type]
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
            
            #check for divByZero
            if(slope == 0):
                ToEs[[k],[i]]=year2
            else:
                #calculate ToE
                Time_of_Emergence = year1 + Standard_Deviation/np.abs(slope)
                if(Time_of_Emergence>year2):
                    ToEs[[k],[i]]=year2
                else:
                    ToEs[[k],[i]] = Time_of_Emergence
            
    lb.title(Precipitation_Type)
    lb.ylabel("mm/day")
    lb.xlabel("Year")
    lb.show()       
    #close file
    ncFile.close()

def generate_WRF_TASMAX_ToE_Data(file, Temperature_Type):
    
    #check input before working on files
    
    if(Temperature_Type !="T2MAX90" and Temperature_Type !="T2MAXx"):
        print("invalid Temperature data type")
        return 1 
    
    #open netcdf file
    
    ncFile = Dataset(file,"r+", format="NETCDF4")
    
        
    # get the data to plot
    TEMP=ncFile.variables[Temperature_Type][:]
    
    #get time axis and grid values
    year=shape(TEMP)[0]
    year0=1970
    year=arange(year)+year0
    year1=1999
    year2=2099
    ny = shape(TEMP)[1] # south-north - 123 - wrf / south-north - 21 - gcm
    nx = shape(TEMP)[2] # west-east - 162 - wrf / west-east - 33 - gcm
    
                
    TEMP=TEMP-273.15 #convert to Kelvin
    lats = ncFile.variables["XLAT"][:]
    lons = ncFile.variables["XLONG"][:]
        
    #Prepare Variables for ToE Calculation
    if(ncFile.variables.__str__().find("regressionValues_Slope_for_" + Temperature_Type)==-1):
        regressionValues_Slopes = ncFile.createVariable("regressionValues_Slope_for_" + Temperature_Type,
                                                        "float32" , ("south_north","west_east",))
        regressionValues_Slopes.units = "kg per m^2 per year"
        regressionValues_Yints =  ncFile.createVariable("regressionValues_Yint_for_" + Temperature_Type,
                                                        "float32", ("south_north","west_east",))
        regressionValues_Yints.units = "kg per m^2"
    if(ncFile.variables.__str__().find("Standard_Deviations_for_" + Temperature_Type)==-1):
        std = ncFile.createVariable("Standard_Deviations_for_" + Temperature_Type, "float32" , ("south_north","west_east",))
        std.units = "kg per m^2"
    if(ncFile.variables.__str__().find("Means_for_" + Temperature_Type)==-1):
        mean = ncFile.createVariable("Means_for_" + Temperature_Type, "float32" , ("south_north","west_east",))
        mean.units = "kg per m^2" 
    #define ToE variable at each grid point
    if(ncFile.variables.__str__().find("ToE_for_" + Temperature_Type)==-1):
        ToE = ncFile.createVariable("ToE_for_" + Temperature_Type, "float32", ("south_north","west_east",))
        ToE.units = "years"
    else:
        print("ToE data already generated for " + Temperature_Type)
        return 1
    #iterate throught the grid points and assign the new variables their respective values
    
    Slopes   = ncFile.variables["regressionValues_Slope_for_" + Temperature_Type]
    Yints    = ncFile.variables["regressionValues_Yint_for_" + Temperature_Type]
    stds     = ncFile.variables["Standard_Deviations_for_" + Temperature_Type]
    means    = ncFile.variables["Means_for_" + Temperature_Type]
    ToEs     = ncFile.variables["ToE_for_" + Temperature_Type]
    
    for k in range(ny):
        for i in range(nx):
            
            #find the mean and standard deviations
            Standard_Deviation=np.std(TEMP[:,[k],[i]][0:31])
            mn=np.mean(TEMP[:,[k],[i]][0:31])
            
            #get mean plus std and record it in file
            means[[k],[i]]=mn
            stds[[k],[i]]=Standard_Deviation
            
            #make two arrays to find a linear regression for a single grid point
            x = year[30:130]
            y = TEMP[:,[k],[i]][30:130]
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
            
            #check for divByZero
            if(slope == 0):
                ToEs[[k],[i]]=year2
            else:
                #calculate ToE
                Time_of_Emergence = year1 + Standard_Deviation/np.abs(slope)
                if(Time_of_Emergence>year2):
                    ToEs[[k],[i]]=year2
                else:
                    ToEs[[k],[i]] = Time_of_Emergence
            
    
    #close file
    ncFile.close()    
    
def generate_WRF_PREC_ToE_Data(file, Precipitation_Type):
    
    #check input before working on files
    
    if(Precipitation_Type !="PREC95" and Precipitation_Type !="PRECx"):
        print("invalid Precipitation data type")
        return 1 
    
    #open netcdf file
    
    ncFile = Dataset(file,"r+", format="NETCDF4")
    PREC=ncFile.variables[Precipitation_Type][:]
    
    #get time axis and grid values
    year=shape(PREC)[0]
    year0=1970
    year=arange(year)+year0
    year1=1999
    year2=2099
    ny = shape(PREC)[1] # south-north - 123
    nx = shape(PREC)[2] # west-east - 162
        
    #check boundary points
    #LOGIC: the values at the boundaries will be rewritten to have values of thier neighbors
    #       towards the inside of the grid.
    #        -------    AABCDEE
    #        -ABCDE-    AABCDEE
    #        -FGHIJ- => FFGHIJJ
    #        -KLMNO-    KKLMNOO
    #        -------    KKLMNOO
    #
    for t in range(year2-year0+1):
        for k in range(ny):
            for i in range(nx):
            #corners  
                if(k==0 and i==0): #top left
                    PREC[[t],[k],[i]]=PREC[[t],[k+1],[i+1]]
                    
                elif(k==ny-1 and i==nx-1): # bottom right
                    PREC[[t],[k],[i]]=PREC[[t],[k-1],[i-1]]
                    
                elif(k==ny-1 and i==0): #bottom left
                    PREC[[t],[k],[i]]=PREC[[t],[k-1],[i+1]]
                    
                elif(k==0 and i==nx-1): #top right
                    PREC[[t],[k],[i]]=PREC[[t],[k+1],[i-1]]          
            #edges      
                elif(k==0 and i!=0 and i!=nx-1): #top edge
                    PREC[[t],[k],[i]]=PREC[[t],[k+1],[i]]
                    
                elif(k!=0 and k!=ny-1 and i==0): #left edge
                    PREC[[t],[k],[i]]=PREC[[t],[k],[i+1]]
                    
                elif(k==ny-1 and i!=0 and i!=nx-1): #bottom edge
                    PREC[[t],[k],[i]]=PREC[[t],[k-1],[i]]
                    
                elif(k!=0 and k!=ny-1 and i==nx-1): #rigth edge
                    PREC[[t],[k],[i]]=PREC[[t],[k],[i-1]]            
                
    PREC=PREC*24*60*60 # mm/day
    
    lats = ncFile.variables["XLAT"][:]
    lons = ncFile.variables["XLONG"][:]
    
    
    #Prepare Variables for ToE Calculation
    if(ncFile.variables.__str__().find("regressionValues_Slope_for_" + Precipitation_Type)==-1):
        regressionValues_Slopes = ncFile.createVariable("regressionValues_Slope_for_" + Precipitation_Type,
                                                        "float32" , ("south_north","west_east",))
        regressionValues_Slopes.units = "mm per day per year"
        regressionValues_Yints =  ncFile.createVariable("regressionValues_Yint_for_" + Precipitation_Type,
                                                        "float32", ("south_north","west_east",))
        regressionValues_Yints.units = "mm/day"
    if(ncFile.variables.__str__().find("Standard_Deviations_for_" + Precipitation_Type)==-1):
        std = ncFile.createVariable("Standard_Deviations_for_" + Precipitation_Type, "float32" , ("south_north","west_east",))
        std.units = "mm/day"
    if(ncFile.variables.__str__().find("Means_for_" + Precipitation_Type)==-1):
        mean = ncFile.createVariable("Means_for_" + Precipitation_Type, "float32" , ("south_north","west_east",))
        mean.units = "mm/day" 
    #define ToE variable at each grid point
    if(ncFile.variables.__str__().find("ToE_for_" + Precipitation_Type)==-1):
        ToE = ncFile.createVariable("ToE_for_" + Precipitation_Type, "float32", ("south_north","west_east",))
        ToE.units = "years"
    else:
        print("ToE data already generated for " + Precipitation_Type)
        return 1
    #iterate throught the grid points and assign the new variables their respective values
    
    Slopes   = ncFile.variables["regressionValues_Slope_for_" + Precipitation_Type]
    Yints    = ncFile.variables["regressionValues_Yint_for_" + Precipitation_Type]
    stds     = ncFile.variables["Standard_Deviations_for_" + Precipitation_Type]
    means    = ncFile.variables["Means_for_" + Precipitation_Type]
    ToEs     = ncFile.variables["ToE_for_" + Precipitation_Type]
    
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
            #check for divByZero
            if(slope == 0):
                ToEs[[k],[i]]=year2
            else:
                #calculate ToE
                Time_of_Emergence = year1 + Standard_Deviation/np.abs(slope)
                if(Time_of_Emergence>year2):
                    ToEs[[k],[i]] = year2
                else:
                    ToEs[[k],[i]] = Time_of_Emergence
            
    #close file
    ncFile.close()
    
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
def interpolate_ToE(file,grid_File, data_Type):        
            
    data=Dataset(file, "r+", format="NETCDF4")
    toe=data.variables["ToE_for_" + data_Type][:] 
    lon=data.variables["lon"][:] - 360
    lat=data.variables["lat"][:] 
    
    # Here use x as lon and y as lat and Z as data field from GCM
    x = lon
    y = lat
    Z = toe
    # here X2 and Y2 (upper case) are lon and lat from WRF grid.nc
    # note that lat and lon are 2d arrays 
    #x2 = linspace(0, 1, 123)
    #y2 = linspace(0, 1, 162)
    #X2,Y2 = meshgrid(x2,y2) # from grid.nc
    
    wdata=Dataset(grid_File, "r", format="NETCDF4")
    wrf_lon=wdata.variables["XLONG"][:] 
    wrf_lat=wdata.variables["XLAT"][:] 
    if(data.dimensions.__str__().find("wrf-latitude")==-1):
        data.createDimension( "wrf-latitude" , size=162)
        data.createDimension( "wrf-longitude" , size=123)
    interpolated_ToE = data.createVariable("Interpolated ToE data based on " + data_Type, "float32" , ("wrf-longitude","wrf-latitude",))              
    interpolated_ToE.units = "years" 
    ToEs = data.variables["Interpolated ToE data based on " + data_Type]
    wrf_interp = interp2d(lon, lat, toe, kind='cubic')

    # Since X2 and Y2 are 2-d arrays for irregular grid need to do point by point
    #    **THERE MAY BE A BETTER WAY?

    wrf_pr=zeros_like(wrf_lon)
    for i in arange(wrf_lon.shape[0]):
        for j in arange(wrf_lon.shape[1]):
            wrf_pr[i,j]=wrf_interp(wrf_lon[i,j],wrf_lat[i,j])
            ToEs[[i],[j]]=wrf_pr[i,j]
        

    # plot coarse and fine grids
    # LON, LAT = meshgrid(lon,lat)

    # fig, ax = plt.subplots(nrows=1, ncols=2)
    # ax[0].pcolormesh(LON,LAT,toe, shading="auto")

    # ax[1].pcolormesh(wrf_lon,wrf_lat,wrf_pr, shading="auto")
    # head, tail = os.path.split(file)
    # x=tail.split("_")
    # plt.figure(figsize=(15,10))
    # WRFplot(wrf_pr,wrf_lat,wrf_lon, amin(wrf_pr),amax(wrf_pr), x[0] + " " + x[1] + " Based On " + data_Type, "ToE in years" , "RdYlBu_r")
    # plt.show()
    # plt.show()

    #close files         
    wdata.close()
    data.close()         
        
def get_Medain_ToE(file, data_Type):
    data=Dataset(file, "r", format="NETCDF4")
    toe=data.variables["ToE_for_" + data_Type][:] 
    head, tail = os.path.split(file)
    x=tail.split("_")
    
    if(data_Type !="PREC95" and data_Type !="PRECx" and
       data_Type !="T2MAX90" and data_Type !="T2MAXx" and 
       data_Type !="prx" and data_Type !="pr95" and 
       data_Type !="tasmaxx" and data_Type != "tasmax90"):
        print("invalid data type")
        return 1 
    #turn 2d masked array into 1 d array so numpy can sort and find the median
    array = toe.flatten()
    print("The median ToE given by " + x[0] + " " + x[1] + " Based On " + data_Type + " is:", np.median(array))
    
    data.close()    
    return ((x[0] + " " + x[1]) , (data_Type) , (np.median(array)))
#execute
import os



#for each wrf and gcm file, generate toe data for each data type and get interpolated gcm data
# (wrf, gcm) (data_Type)
gcm_t_Data_Types = ["tasmax90", "tasmaxx"]
gcm_p_Data_Types = ["pr95", "prx"]
wrf_t_Data_Types = ["T2MAX90", "T2MAXx"]
wrf_p_Data_Types = ["PREC95", "PRECx"]
wrf_Folder = r"Netcdf_Files" + "\\" +"wrf_Netcdf_Files"
gcm_Folder = r"Netcdf_Files" + "\\" +"gcm_Netcdf_Files"
#for each wrf and gcm file, generate toe data for each data type
for wrf_File in os.listdir(wrf_Folder): 
    wrf_head, wrf_tail = os.path.split(wrf_File) #tail gives the file name and type
    if(wrf_tail.__str__().find("T2MAX_extr")>=0):
        for temperature_Type in wrf_t_Data_Types:
            generate_WRF_TASMAX_ToE_Data("Netcdf_Files" + "\\" + "wrf_Netcdf_Files" + "\\" + wrf_tail, temperature_Type)
    if(wrf_tail.__str__().find("PREC_extr")>=0):
        for precipitation_Type in wrf_p_Data_Types:
            generate_WRF_PREC_ToE_Data("Netcdf_Files" + "\\" + "wrf_Netcdf_Files" + "\\" + wrf_tail, precipitation_Type)
            
for gcm_File in os.listdir(gcm_Folder):
    gcm_head, gcm_tail = os.path.split(gcm_File) 
    if(gcm_tail.__str__().find("tasmax_extr")>=0):
        for temperature_Type in gcm_t_Data_Types:
            generate_GCM_TASMAX_ToE_Data("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + gcm_tail, temperature_Type)
            interpolate_ToE("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + gcm_tail,"wrf_grid.nc", temperature_Type)
    if(gcm_tail.__str__().find("pr_extr")>=0):
        for precipitation_Type in gcm_p_Data_Types:
            generate_GCM_PREC_ToE_Data("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + gcm_tail, precipitation_Type)
            interpolate_ToE("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + gcm_tail,"wrf_grid.nc", precipitation_Type)


#for each data type and model calculate and store the median time of emergence.
median_File = open("Median_Toe.dat", "w")
median_File.write("Model Name               data type        Median ToE")  # write table header
Models = []
data_Types = []
median_ToEs = []
for wrf_File in os.listdir(wrf_Folder): 
    wrf_head, wrf_tail = os.path.split(wrf_File) #tail gives the file name and type
    if(wrf_tail.__str__().find("T2MAXextr")>=0):
        for temperature_Type in wrf_t_Data_Types:
            median = get_Medain_ToE("Netcdf_Files" + "\\" + "wrf_Netcdf_Files" + "\\" + wrf_tail, temperature_Type)
            Models.append(median[0])
            data_Types.append(median[1])
            median_ToEs.append(median[2])
    if(wrf_tail.__str__().find("PRECextr")>=0):
        for precipitation_Type in wrf_p_Data_Types:
            median = get_Medain_ToE("Netcdf_Files" + "\\" + "wrf_Netcdf_Files" + "\\" + wrf_tail, precipitation_Type)
            Models.append(median[0])
            data_Types.append(median[1])
            median_ToEs.append(median[2])
for gcm_File in os.listdir(gcm_Folder):
    gcm_head, gcm_tail = os.path.split(gcm_File) 
    if(gcm_tail.__str__().find("tasmaxextr")>=0):
        for temperature_Type in gcm_t_Data_Types:
            median = get_Medain_ToE("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + gcm_tail, temperature_Type)
            Models.append(median[0])
            data_Types.append(median[1])
            median_ToEs.append(median[2])
    if(gcm_tail.__str__().find("prextr")>=0):
        for precipitation_Type in gcm_p_Data_Types:
            median = get_Medain_ToE("Netcdf_Files" + "\\" + "gcm_Netcdf_Files" + "\\" + gcm_tail, precipitation_Type)
            Models.append(median[0])
            data_Types.append(median[1])
            median_ToEs.append(median[2])

for model, data_Type, median_ToE  in zip(Models, data_Types, median_ToEs):
    median_File.write("\n")
    median_File.write(model + "    " + data_Type + "        " + str(median_ToE))
 
