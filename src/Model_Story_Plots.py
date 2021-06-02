# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 16:55:36 2021

@author: Kedar
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 09:49:33 2021

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

from PIL import Image
import matplotlib.image as mpimg
import sys

def Generate_Story_Plots():
    models = ["access1.3","bcc","canesm2","ccsm4","csiro","fgoals","gfdl","giss","miroc5","mri","noresm1"]
    Ensemble_Slope_Folder = r"Ensemble Slope Plots" 
    Ensemble_SDV_Folder = r"Ensemble SDV Plots"
    Ensemble_TOE_Folder = r"Ensemble TOE Plots"
    
    
    t_Data_Types = ["tasmax90"]
    p_Data_Types = ["pr95"]
    
    
    for model in models:
        for temperature_Type in t_Data_Types:
             Images = []
             wrf_Data_Type = " "
             if(temperature_Type == "tasmax90"):
                 wrf_Data_Type = "T2MAX90"
             elif(temperature_Type == "tasmaxx"):
                 wrf_Data_Type = "T2MAXx"
             elif(temperature_Type !="tasmax90" and temperature_Type != "tasmaxx"):
                 print("invalid temperature_Type")
                 return 1
             for TOE_IMG in os.listdir(Ensemble_TOE_Folder):
                 TOE_IMG_HEAD , TOE_IMG_TAIL = os.path.split(TOE_IMG)
                 if(TOE_IMG_TAIL.__str__().find(model)>=0 and TOE_IMG_TAIL.__str__().find(temperature_Type)>=0):
                     Images.append(Image.open(Ensemble_TOE_Folder + "\\"+ TOE_IMG))
                     
             for SDV_IMG in os.listdir(Ensemble_SDV_Folder):
                 SDV_IMG_HEAD , SDV_IMG_TAIL = os.path.split(SDV_IMG)
                 if(SDV_IMG_TAIL.__str__().find(model)>=0 and SDV_IMG_TAIL.__str__().find(temperature_Type)>=0):
                     Images.append(Image.open(Ensemble_SDV_Folder + "\\"+ SDV_IMG))
                     
             for SLOPE_IMG in os.listdir(Ensemble_Slope_Folder):
                 SLOPE_IMG_HEAD , SLOPE_IMG_TAIL = os.path.split(SLOPE_IMG)
                 if(SLOPE_IMG_TAIL.__str__().find(model)>=0 and SLOPE_IMG_TAIL.__str__().find(temperature_Type)>=0):
                     Images.append(Image.open(Ensemble_Slope_Folder + "\\"+ SLOPE_IMG))
                     
         
             
             widths, heights = zip(*(i.size for i in Images))
            
             total_height = sum(heights)
             max_width = max(widths)
            
             new_im = Image.new('RGBA', (max_width, total_height))
            
             
             y_offset = 0
             for im in Images:
                 new_im.paste(im, (0,y_offset))
                 y_offset += im.size[1]
            
             new_im.save("Model Story Plots"+ "\\" + model + " based on " + temperature_Type+".png")
             
        for precipitation_Type in p_Data_Types:
             Images = []
             wrf_Data_Type = " "
             if(precipitation_Type == "pr95"):
                 wrf_Data_Type = "PREC95"
             elif(precipitation_Type == "prx"):
                 wrf_Data_Type = "PRECx"
             elif(precipitation_Type != "prx" and precipitation_Type != "pr95"):
                 print("invalid precipitation_Type")
                 return 1
             for TOE_IMG in os.listdir(Ensemble_TOE_Folder):
                 TOE_IMG_HEAD , TOE_IMG_TAIL = os.path.split(TOE_IMG)
                 if(TOE_IMG_TAIL.__str__().find(model)>=0 and TOE_IMG_TAIL.__str__().find(precipitation_Type)>=0):
                     Images.append(Image.open(Ensemble_TOE_Folder + "\\"+ TOE_IMG))
                     
             for SDV_IMG in os.listdir(Ensemble_SDV_Folder):
                 SDV_IMG_HEAD , SDV_IMG_TAIL = os.path.split(SDV_IMG)
                 if(SDV_IMG_TAIL.__str__().find(model)>=0 and SDV_IMG_TAIL.__str__().find(precipitation_Type)>=0):
                     Images.append(Image.open(Ensemble_SDV_Folder + "\\"+ SDV_IMG))
                     
             for SLOPE_IMG in os.listdir(Ensemble_Slope_Folder):
                 SLOPE_IMG_HEAD , SLOPE_IMG_TAIL = os.path.split(SLOPE_IMG)
                 if(SLOPE_IMG_TAIL.__str__().find(model)>=0 and SLOPE_IMG_TAIL.__str__().find(precipitation_Type)>=0):
                     Images.append(Image.open(Ensemble_Slope_Folder + "\\"+ SLOPE_IMG))
        
     
             widths, heights = zip(*(i.size for i in Images))
            
             total_height = sum(heights)
             max_width = max(widths)
            
             new_im = Image.new('RGBA', (max_width, total_height))
            
             y_offset = 0
             for im in Images:
                 new_im.paste(im, (0,y_offset))
                 y_offset += im.size[1]

             new_im.save("Model Story Plots"+ "\\" + model + " based on " + precipitation_Type+".png")
             

#execute
Generate_Story_Plots()