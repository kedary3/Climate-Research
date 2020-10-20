# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:35:10 2020

@author: Kedar Yadav
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

def generate_WRF_TASMAX_ToE_Data(file, Temperature_Type):
    
    #check input before working on files
    
    if(Temperature_Type !="T2MAX90" and Temperature_Type !="T2MAXx"):
        print("invalid Temperature data type")
        return 1 
    
    #open netcdf file
    
    ncFile = Dataset(file,"r+", format="NETCDF4")
    
        
    # get the data to plot
    TEMP=ncFile.variables[Temperature_Type][:]-273.15 #convert to Kelvin
    lats = ncFile.variables["XLAT"][:]
    lons = ncFile.variables["XLONG"][:]
    
    #get time axis and grid values
    year=shape(TEMP)[0]
    year0=1970
    year=arange(year)+year0
    year1=1999
    year2=2099
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
    lb.figure()
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
            #Plot the linear regression
            lb.plot(x_lin_reg, y_lin_reg, c = 'black')
            #plot the data at that grid point over time
            lb.plot(year,TEMP[:,[k],[i]])
            #check for divByZero
            if(slope == 0):
                ToEs[[k],[i]]=3000
            else:
                #calculate ToE
                if(slope<0):
                    Time_of_Emergence = (Standard_Deviation-mn-Yint)/slope
                    if(Time_of_Emergence>3000):
                        ToEs[[k],[i]] = 3000
                    else:
                        ToEs[[k],[i]] = Time_of_Emergence
                else:
                    Time_of_Emergence = (Standard_Deviation+mn-Yint)/slope
                    if(Time_of_Emergence>3000):
                        ToEs[[k],[i]] = 3000
                    else:
                        ToEs[[k],[i]] = Time_of_Emergence
            
    lb.title(Temperature_Type)
    lb.ylabel("Tmax (deg-C)")
    lb.xlabel("Year")
    lb.show()       
    #close file
    ncFile.close()
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
            #Plot the linear regression
            lb.plot(x_lin_reg, y_lin_reg, c = 'black')
            #plot the data at that grid point over time
            lb.plot(year,TEMP[:,[k],[i]])
            #check for divByZero
            if(slope == 0):
                ToEs[[k],[i]]=3000
            else:
                #calculate ToE
                if(slope<0):
                    Time_of_Emergence = (Standard_Deviation-mn-Yint)/slope
                    if(Time_of_Emergence>3000):
                        ToEs[[k],[i]] = 3000
                    else:
                        ToEs[[k],[i]] = Time_of_Emergence
                else:
                    Time_of_Emergence = (Standard_Deviation+mn-Yint)/slope
                    if(Time_of_Emergence>3000):
                        ToEs[[k],[i]] = 3000
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
    
    if(Precipitation_Type !="prx" and Precipitation_Type !="pr95" and Precipitation_Type !="prn"):
        print("invalid Precipitation data type")
        return 1 
    
    #open netcdf file
    
    ncFile = Dataset(file,"r+", format="NETCDF4")
    
        
    # get the data to plot
    PREC=ncFile.variables[Precipitation_Type][:]
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
    #splice the variable array along the time axis, separating model and historical data
    #31 is chosen because python is exclusive on the second range number
    Historical_PREC = PREC[:,0,0][0:31]
    plt.plot(year[0:31],Historical_PREC)
    Model_PREC = PREC[:,0,0][30:130]
    plt.plot(year[30:130],Model_PREC)
    plt.show()
    
    #Prepare Variables for ToE Calculation
    if(ncFile.variables.__str__().find("regressionValues_Slope_for_" + Precipitation_Type)==-1):
        regressionValues_Slopes = ncFile.createVariable("regressionValues_Slope_for_" + Precipitation_Type,
                                                        "float32" , ("lat","lon",))
        regressionValues_Slopes.units = "kg per m^2 per year"
        regressionValues_Yints =  ncFile.createVariable("regressionValues_Yint_for_" + Precipitation_Type,
                                                        "float32", ("lat","lon",))
        regressionValues_Yints.units = "kg per m^2"
    if(ncFile.variables.__str__().find("Standard_Deviations_for_" + Precipitation_Type)==-1):
        std = ncFile.createVariable("Standard_Deviations_for_" + Precipitation_Type, "float32" , ("lat","lon",))
        std.units = "kg per m^2"
    if(ncFile.variables.__str__().find("Means_for_" + Precipitation_Type)==-1):
        mean = ncFile.createVariable("Means_for_" + Precipitation_Type, "float32" , ("lat","lon",))
        mean.units = "kg per m^2" 
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
            #Plot the linear regression
            lb.plot(x_lin_reg, y_lin_reg, c = 'black')
            #plot the data at that grid point over time
            lb.plot(year,PREC[:,[k],[i]])
            #check for divByZero
            if(slope == 0):
                ToEs[[k],[i]]=3000
            else:
                #calculate ToE
                if(slope<0):
                    Time_of_Emergence = (Standard_Deviation-mn-Yint)/slope
                    if(Time_of_Emergence>3000):
                        ToEs[[k],[i]] = 3000
                    else:
                        ToEs[[k],[i]] = Time_of_Emergence
                else:
                    Time_of_Emergence = (Standard_Deviation+mn-Yint)/slope
                    if(Time_of_Emergence>3000):
                        ToEs[[k],[i]] = 3000
                    else:
                        ToEs[[k],[i]] = Time_of_Emergence
            
    lb.title(Precipitation_Type)
    lb.ylabel("kg/s*m^2")
    lb.xlabel("Year")
    lb.show()       
    #close file
    ncFile.close()
    
def generate_WRF_PREC_ToE_Data(file, Precipitation_Type):
    
    #check input before working on files
    
    if(Precipitation_Type !="PREC95" and Precipitation_Type !="PRECx"):
        print("invalid Precipitation data type")
        return 1 
    
    #open netcdf file
    
    ncFile = Dataset(file,"r+", format="NETCDF4")
    
        
    # get the data to plot
    PREC=ncFile.variables[Precipitation_Type][:]
    lats = ncFile.variables["XLAT"][:]
    lons = ncFile.variables["XLONG"][:]
    
    #get time axis and grid values
    year=shape(PREC)[0]
    year0=1970
    year=arange(year)+year0
    year1=1999
    year2=2099
    ny = shape(PREC)[1] # south-north - 123
    nx = shape(PREC)[2] # west-east - 162
    #splice the variable array along the time axis, separating model and historical data
    #31 is chosen because python is exclusive on the second range number
    Historical_PREC = PREC[:,0,0][0:31]
    plt.plot(year[0:31],Historical_PREC)
    Model_PREC = PREC[:,0,0][30:130]
    plt.plot(year[30:130],Model_PREC)
    plt.show()
    
    #Prepare Variables for ToE Calculation
    if(ncFile.variables.__str__().find("regressionValues_Slope_for_" + Precipitation_Type)==-1):
        regressionValues_Slopes = ncFile.createVariable("regressionValues_Slope_for_" + Precipitation_Type,
                                                        "float32" , ("south_north","west_east",))
        regressionValues_Slopes.units = "kg per m^2 per year"
        regressionValues_Yints =  ncFile.createVariable("regressionValues_Yint_for_" + Precipitation_Type,
                                                        "float32", ("south_north","west_east",))
        regressionValues_Yints.units = "kg per m^2"
    if(ncFile.variables.__str__().find("Standard_Deviations_for_" + Precipitation_Type)==-1):
        std = ncFile.createVariable("Standard_Deviations_for_" + Precipitation_Type, "float32" , ("south_north","west_east",))
        std.units = "kg per m^2"
    if(ncFile.variables.__str__().find("Means_for_" + Precipitation_Type)==-1):
        mean = ncFile.createVariable("Means_for_" + Precipitation_Type, "float32" , ("south_north","west_east",))
        mean.units = "kg per m^2" 
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
            #check for divByZero
            if(slope == 0):
                ToEs[[k],[i]]=3000
            else:
                #calculate ToE
                if(slope<0):
                    Time_of_Emergence = (Standard_Deviation-mn-Yint)/slope
                    if(Time_of_Emergence>3000):
                        ToEs[[k],[i]] = 3000
                    else:
                        ToEs[[k],[i]] = Time_of_Emergence
                else:
                    Time_of_Emergence = (Standard_Deviation+mn-Yint)/slope
                    if(Time_of_Emergence>3000):
                        ToEs[[k],[i]] = 3000
                    else:
                        ToEs[[k],[i]] = Time_of_Emergence
            
    lb.title(Precipitation_Type)
    lb.ylabel("kg/s*m^2")
    lb.xlabel("Year")
    lb.show()       
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
    toe=data.variables["ToE_for_" + data_Type][:] #* 24*60*60 # mm/day
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
        LON, LAT = meshgrid(lon,lat)
    
        fig, ax = plt.subplots(nrows=1, ncols=2)
        ax[0].pcolormesh(LON,LAT,toe, shading="auto")
    
        ax[1].pcolormesh(wrf_lon,wrf_lat,wrf_pr, shading="auto")
        head, tail = os.path.split(file)
        x=tail.split("_")
        plt.figure(figsize=(15,10))
        WRFplot(wrf_pr,wrf_lat,wrf_lon, amin(wrf_pr),amax(wrf_pr), x[0] + " " + x[1] + " Based On " + data_Type, "ToE in years" , "RdYlBu_r")
        plt.show()
    
    
        plt.show()
    
    #close files         
    wdata.close()
    data.close()         
        
#execute
import os
folder = r'Netcdf_Files'
#for each wrf and gcm file, generate toe data for each data type
for file in os.scandir(folder):
    head, tail = os.path.split(file) #tail gives the file name and type
    if(file.__str__().find("prextr")!=-1): #if file has precipitation data
        
        generate_GCM_PREC_ToE_Data("Netcdf_Files" + "\\" + tail, "prx")
        interpolate_ToE("Netcdf_Files" + "\\" + tail,"wrf_grid.nc", "prx")
        
        generate_GCM_PREC_ToE_Data("Netcdf_Files" + "\\" + tail, "pr95")
        interpolate_ToE("Netcdf_Files" + "\\" + tail,"wrf_grid.nc", "pr95")
        
        generate_GCM_PREC_ToE_Data("Netcdf_Files" + "\\" + tail, "prn")
        interpolate_ToE("Netcdf_Files" + "\\" + tail,"wrf_grid.nc", "prn")
        
    if(file.__str__().find("tasmaxextr")!=-1): #if file has tempature data
        
        generate_GCM_TASMAX_ToE_Data("Netcdf_Files" + "\\" + tail, "tasmax90")
        interpolate_ToE("Netcdf_Files" + "\\" + tail,"wrf_grid.nc", "tasmax90")
        
        generate_GCM_TASMAX_ToE_Data("Netcdf_Files" + "\\" + tail, "tasmaxx")
        interpolate_ToE("Netcdf_Files" + "\\" + tail,"wrf_grid.nc", "tasmaxx")
        
        generate_GCM_TASMAX_ToE_Data("Netcdf_Files" + "\\" + tail, "tasmaxn")
        interpolate_ToE("Netcdf_Files" + "\\" + tail,"wrf_grid.nc", "tasmaxn")
    if(file.__str__().find("T2MAXextr")!=-1): #if file has tempature data
    
        generate_WRF_TASMAX_ToE_Data("Netcdf_Files" + "\\" + tail, "T2MAX90")
        generate_WRF_TASMAX_ToE_Data("Netcdf_Files" + "\\" + tail, "T2MAXx")
        
    if(file.__str__().find("PRECextr")!=-1): #if file has precipitation data
    
        generate_WRF_PREC_ToE_Data("Netcdf_Files" + "\\" + tail, "PREC95")
        generate_WRF_PREC_ToE_Data("Netcdf_Files" + "\\" + tail, "PRECx")