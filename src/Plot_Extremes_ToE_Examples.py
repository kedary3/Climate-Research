# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 20:07:39 2020

@author: Kedar Yadav
"""
# #to dos
#     # plot of two extremes to demonstrate two different toEs and add labels
#     # time of emergence nc create std and mn variables
#     # add variables to model netcdf
      # add metadata to variable names of time of emergence
from netCDF4 import Dataset
from numpy import (
        linspace,array, log,exp,sin,cos,sqrt, pi,e, 
        zeros, ones, amin,amax, argmax, arange
        )
from matplotlib.pyplot import plot, show
import numpy as np
import matplotlib.pyplot as plt
# function to generate an extremes graph
def main():
    #retrieve data
    ds1= Dataset("ccsm4_1970-2006_TX90th.nc",
    "r", format="NETCDF4")
    ds2= Dataset("ccsm4_2006-2099_TX90th.nc",
    "r",format="NETCDF4")
    t90h=ds1.variables["T90"][:] - 273.15 # convert to deg-C
    t90m=ds2.variables["T90"][:] - 273.15 
    mn=ds1.variables["Means"][:]
    std=ds1.variables["Standard_Deviations"][:]
    slope=ds2.variables["regressionValues_Slope"][:]
    intercept=ds2.variables["regressionValues_Yint"][:]
    ToE=ds2.variables["ToE"][:]
    ntimeh = np.shape(t90h)[0] # year index
    ny = np.shape(t90h)[1] # south-north - 123
    nx = np.shape(t90h)[2] # west-east - 162
    ntimem= np.shape(t90m)[0]
    year0=1970
    year1=2006
    yearh = arange(ntimeh)+year0
    yearm = arange(ntimem)+year1
    plt.title("90th Percentile Tmax Vs. Year: Demonstrating Times of Emergence")
    plt.ylabel("Tmax (deg-C)")
    plt.xlabel("Year")
    #get two extreme ToEs
    low=findLowToE("ccsm4_2006-2099_TX90th.nc")
    high=findHighToE("ccsm4_2006-2099_TX90th.nc")
    #plot each extreme
    plt.plot(yearh,t90h[:,[low[0]],[low[1]]]) 
    plt.plot(yearm,t90m[:,[low[0]],[low[1]]])
    plt.plot(yearh,t90h[:,[high[0]],[high[1]]]) 
    plt.plot(yearm,t90m[:,[high[0]],[high[1]]])
    #Plot the mean and standard deviation lines for low ToE
    x1= [year0,year1+ntimem] #extra year distance for readability
    y1= [mn[[low[0]],[low[1]]]+std[[low[0]],[low[1]]],mn[[low[0]],[low[1]]]+std[[low[0]],[low[1]]]]
    y2= [mn[[low[0]],[low[1]]]-std[[low[0]],[low[1]]],mn[[low[0]],[low[1]]]-std[[low[0]],[low[1]]]]
    y3= [mn[[low[0]],[low[1]]],mn[[low[0]],[low[1]]]]
    plt.plot(x1,y1,marker="o")
    plt.plot(x1,y2,marker="o")
    plt.plot(x1,y3,marker="o")
    #Plot the mean and standard deviation lines for high ToE
    x1= [year0,year1+ntimem] #extra year distance for readability
    y1= [mn[[high[0]],[high[1]]]+std[[high[0]],[high[1]]],mn[[high[0]],[high[1]]]+std[[high[0]],[high[1]]]]
    y2= [mn[[high[0]],[high[1]]]-std[[high[0]],[high[1]]],mn[[high[0]],[high[1]]]-std[[high[0]],[high[1]]]]
    y3= [mn[[high[0]],[high[1]]],mn[[high[0]],[high[1]]]]
    plt.plot(x1,y1,marker="o")
    plt.plot(x1,y2,marker="o")
    plt.plot(x1,y3,marker="o")
    #Plot the regression line for the low ToE
    x_lin_reg = range(year1, year1+ntimem)
    plt.plot(x_lin_reg, slope[[low[0]],[low[1]]]*x_lin_reg + intercept[[low[0]],[low[1]]], c= "black")
    #Plot the regression line for the high ToE
    plt.plot(x_lin_reg, slope[[high[0]],[high[1]]]*x_lin_reg + intercept[[high[0]],[high[1]]], c= "black")
    #Plot the times of emergence for the low
    t1,y1 = [ToE[[low[0]],[low[1]]],ToE[[low[0]],[low[1]]]],[mn[[low[0]],[low[1]]]-std[[high[0]],[low[1]]],mn[[low[0]],[low[1]]]+std[[low[0]],[low[1]]]]
    plt.plot(t1,y1, marker="_", c="black")
    #Plot the times of emergence for the high
    t1,y1 = [ToE[[high[0]],[high[1]]],ToE[[high[0]],[high[1]]]],[mn[[high[0]],[high[1]]]-std[[high[0]],[high[1]]],mn[[high[0]],[high[1]]]+std[[high[0]],[high[1]]]]
    plt.plot(t1,y1, marker="_", c="black")
    #annotate the ToEs
    plt.annotate("Low Time of Emergence", (ToE[[low[0]],[low[1]]],(mn[[low[0]],[low[1]]]-std[[low[0]],[low[1]]]-1)))
    plt.annotate("High Time of Emergence", (ToE[[high[0]],[high[1]]]-50,(mn[[high[0]],[high[1]]]+2*std[[high[0]],[high[1]]]+1)))

    plt.savefig("ToE Extremes")
    plt.show()
    
# find the grid point that has the highest ToE    
def findHighToE(fl):
    ds=Dataset(fl,"r", format="NETCDF4")
    Max_Location=(0,0)
    ToE=ds.variables["ToE"][:]
    ny = np.shape(ToE)[0] # south-north - 123
    nx = np.shape(ToE)[1] # west-east - 162
    ToE_Max=0
    for k in range(ny):
        for i in range(nx):
            if(ToE[[k],[i]]>=ToE_Max):  #check if the new grid point has a higher ToE
                ToE_Max=ToE[[k],[i]]
                Max_Location=(k,i)
            else:
                continue                #if not we want to continue looking through the gridpoints
    ds.close()
    return Max_Location
def findLowToE(fl):
    ds=Dataset(fl,"r", format="NETCDF4")
    Low_Location=(0,0)
    ToE=ds.variables["ToE"][:]
    ny = np.shape(ToE)[0] # south-north - 123
    nx = np.shape(ToE)[1] # west-east - 162
    ToE_Low=3000
    for k in range(ny):
        for i in range(nx):
            if(ToE[[k],[i]]<=ToE_Low):  #check if the new grid point has a lower ToE
                ToE_Low=ToE[[k],[i]]
                Low_Location=(k,i)
            else:
                continue                #if not we want to continue looking through the gridpoints
    ds.close()
    return Low_Location
#execute 
main()
#to does
# move labels to outside graph
# look at new data to try dynamic code
# 