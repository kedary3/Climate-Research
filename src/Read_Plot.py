#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 12:59:55 2020

@author: Kedar Yadav
"""
from Global_Vars import mean as mn
from Global_Vars import std 
import numpy as np
import matplotlib.pyplot as plt 
#from sklearn.metrics import r2_score

# repeat for each file in list
# from netCDF4 import Dataset
# Net_data = Dataset("ccsm4_2006-2099_TX90th.nc", "r" ,format="NETCDF4")
# print(Net_data)
# for name in Net_data.ncattrs():
#      print("Global attr {} = {}".format(name, getattr(Net_data, name)))
                                                      
# Net_data.close
for file in [
        "ccsm4_tx90p_ksea_historicalNat.dat",
        "ccsm4_tx90p_ksea_rcp85.dat"
        ]:

    # Read data file  to big array
    data = np.loadtxt(file)
    
    # First column is year
    year = data[:,0]
    
    # second is percent of days exceeding historic 90th percentile (tx90p)
    tx90p = data[:,1]
    
    # plot the data
    plt.plot(year,tx90p)
    #to dos 
    
    #if the current file is of historical data, find the standard deviation
    if(file.__str__().find("historical")!=-1): 
        
        Standard_Deviation=np.std(tx90p)
        #set globals so they can be accessed
        std=Standard_Deviation 
        mn=np.mean(tx90p)
        #Plot the mean and standard deviation lines
        x1= [year[0],int(year[len(year)-1]+len(year))]
        y1= [mn+(Standard_Deviation),mn+(Standard_Deviation)] 
        y2= [mn-(Standard_Deviation),mn-(Standard_Deviation)]
        y3= [np.mean(tx90p),np.mean(tx90p)]
        plt.plot(x1,y1,marker="o")
        plt.plot(x1,y2,marker="o")
        plt.plot(x1,y3,marker="o")
    #if the current file is a model, find its linear regression and the ToE
    if(file.__str__().find("historical")==-1): 
        x=year
        y=tx90p
        Linear_Regression= np.polyfit(x,y,1)
        predict=np.poly1d(Linear_Regression)
        #slope of regression is lin_reg[0], y intercept is lin_reg[1]
        #R2+score is the accuracy of the linear fit
        #r2_score=r2_score(y, predict(x))
        x_lin_reg = range(int(year[0]), int(year[len(year)-1]))
        y_lin_reg = predict(x_lin_reg)
        #Plot the linear regression
        plt.plot(x_lin_reg, y_lin_reg, c = 'black')
        #Calculate the ToE
        if(Linear_Regression[0]==0):
            Time_of_Emergence="Never"
        else:
            Time_of_Emergence = (mn+(std)-Linear_Regression[1])/Linear_Regression[0]
            print("The Time of Emergence is " + Time_of_Emergence.__str__())
            t1,y1 = [Time_of_Emergence,Time_of_Emergence], [0,np.std(tx90p)]
            #Plot the ToE line
            plt.plot(t1,y1, marker="_")
   
    
    