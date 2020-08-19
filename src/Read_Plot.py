#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 12:59:55 2020

@author: esalathe
"""
import numpy as np
import matplotlib.pyplot as plt 
# repeat for each file in list
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
    Standard_Deviation=np.std(tx90p)
    print(Standard_Deviation)
    # get the standard deviation of historical 
    # alpha= sqrt((1/N)*sum (i=0 to N)(x(i)-mean)^2)
    # linear fit the rcp predictied data
    # calculate the ToE

    
    