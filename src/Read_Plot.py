#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 12:59:55 2020

@author: esalathe
"""
import numpy as np
print(np.version)
# repeat for each file in list
for file in [
        "ccsm4_tx90p_ksea_historicalNat.dat",
        "ccsm4_tx90p_ksea_rcp85.dat"
        ]:

    # Read data file  to big array
    data = numpy.loadtxt(file)
    
    # First column is year
    year = data[:,0]
    
    # second is percent of days exceeding historic 90th percentile (tx90p)
    tx90p = data[:,1]
    
    # plot the data
    numpy.plot(year,tx90p)
    
    