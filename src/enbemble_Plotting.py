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

def Generate_Three_Pane_Slope_Plots():
    models = ["access1.3","bcc","canesm2","ccsm4","csiro","fgoals","gfdl","giss","miroc5","mri","noresm1"]
    WRF_Slope_Folder = r"WRF Slope Plots" 
    GCM_Slope_Folder = r"GCM Slope Plots"
    Delta_Slope_Folder = r"Slope Deltas"
    
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
             for WRF_IMG in os.listdir(WRF_Slope_Folder):
                 WRF_IMG_HEAD , WRF_IMG_TAIL = os.path.split(WRF_IMG)
                 if(WRF_IMG_TAIL.__str__().find(model)>=0 and WRF_IMG_TAIL.__str__().find(wrf_Data_Type)>=0):
                     Images.append(Image.open(WRF_Slope_Folder + "\\"+ WRF_IMG))
             for GCM_IMG in os.listdir(GCM_Slope_Folder): 
                 GCM_IMG_HEAD , GCM_IMG_TAIL = os.path.split(GCM_IMG)
                 if(GCM_IMG_TAIL.__str__().find(model)>=0 and GCM_IMG_TAIL.__str__().find(temperature_Type)>=0):
                     Images.append(Image.open(GCM_Slope_Folder + "\\"+GCM_IMG))
             for Delta_IMG in os.listdir(Delta_Slope_Folder): 
                 Delta_IMG_HEAD , Delta_IMG_TAIL = os.path.split(Delta_IMG)
                 if(Delta_IMG_TAIL.__str__().find(model)>=0 and Delta_IMG_TAIL.__str__().find(temperature_Type)>=0):
                     Images.append(Image.open(Delta_Slope_Folder + "\\"+Delta_IMG))
        
     
             widths, heights = zip(*(i.size for i in Images))
            
             total_width = sum(widths)
             max_height = max(heights)
            
             new_im = Image.new('RGBA', (total_width, max_height))
            
             x_offset = 0
             for im in Images:
                 new_im.paste(im, (x_offset,0))
                 x_offset += im.size[0]

             new_im.save("Ensemble Slope Plots"+ "\\" + "Ensemble Slope Plot with " + model + " based on " + temperature_Type+".png")
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
             for WRF_IMG in os.listdir(WRF_Slope_Folder):
                 WRF_IMG_HEAD , WRF_IMG_TAIL = os.path.split(WRF_IMG)
                 if(WRF_IMG_TAIL.__str__().find(model)>=0 and WRF_IMG_TAIL.__str__().find(wrf_Data_Type)>=0):
                     Images.append(Image.open(WRF_Slope_Folder + "\\"+ WRF_IMG))
             for GCM_IMG in os.listdir(GCM_Slope_Folder): 
                 GCM_IMG_HEAD , GCM_IMG_TAIL = os.path.split(GCM_IMG)
                 if(GCM_IMG_TAIL.__str__().find(model)>=0 and GCM_IMG_TAIL.__str__().find(precipitation_Type)>=0):
                     Images.append(Image.open(GCM_Slope_Folder + "\\"+GCM_IMG))
             for Delta_IMG in os.listdir(Delta_Slope_Folder): 
                 Delta_IMG_HEAD , Delta_IMG_TAIL = os.path.split(Delta_IMG)
                 if(Delta_IMG_TAIL.__str__().find(model)>=0 and Delta_IMG_TAIL.__str__().find(precipitation_Type)>=0):
                     Images.append(Image.open(Delta_Slope_Folder + "\\"+Delta_IMG))
        
     
             widths, heights = zip(*(i.size for i in Images))
            
             total_width = sum(widths)
             max_height = max(heights)
            
             new_im = Image.new('RGBA', (total_width, max_height))
            
             x_offset = 0
             for im in Images:
                 new_im.paste(im, (x_offset,0))
                 x_offset += im.size[0]

             new_im.save("Ensemble Slope Plots"+ "\\" + "Ensemble Slope Plot with "+ model + " based on " + precipitation_Type+".png")


def Generate_Three_Pane_SDV_Plots():
    models = ["access1.3","bcc","canesm2","ccsm4","csiro","fgoals","gfdl","giss","miroc5","mri","noresm1"]
    WRF_SDV_Folder = r"WRF SDV Plots" 
    GCM_SDV_Folder = r"GCM SDV Plots"
    Delta_SDV_Folder = r"SDV Deltas"
    
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
             for WRF_IMG in os.listdir(WRF_SDV_Folder):
                 WRF_IMG_HEAD , WRF_IMG_TAIL = os.path.split(WRF_IMG)
                 if(WRF_IMG_TAIL.__str__().find(model)>=0 and WRF_IMG_TAIL.__str__().find(wrf_Data_Type)>=0):
                     Images.append(Image.open(WRF_SDV_Folder + "\\"+ WRF_IMG))
             for GCM_IMG in os.listdir(GCM_SDV_Folder): 
                 GCM_IMG_HEAD , GCM_IMG_TAIL = os.path.split(GCM_IMG)
                 if(GCM_IMG_TAIL.__str__().find(model)>=0 and GCM_IMG_TAIL.__str__().find(temperature_Type)>=0):
                     Images.append(Image.open(GCM_SDV_Folder + "\\"+GCM_IMG))
             for Delta_IMG in os.listdir(Delta_SDV_Folder): 
                 Delta_IMG_HEAD , Delta_IMG_TAIL = os.path.split(Delta_IMG)
                 if(Delta_IMG_TAIL.__str__().find(model)>=0 and Delta_IMG_TAIL.__str__().find(temperature_Type)>=0):
                     Images.append(Image.open(Delta_SDV_Folder + "\\"+Delta_IMG))
        
     
             widths, heights = zip(*(i.size for i in Images))
            
             total_width = sum(widths)
             max_height = max(heights)
            
             new_im = Image.new('RGBA', (total_width, max_height))
            
             x_offset = 0
             for im in Images:
                 new_im.paste(im, (x_offset,0))
                 x_offset += im.size[0]

             new_im.save("Ensemble SDV Plots"+ "\\" + "Ensemble SDV Plot with " + model + " based on " + temperature_Type+".png")
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
             for WRF_IMG in os.listdir(WRF_SDV_Folder):
                 WRF_IMG_HEAD , WRF_IMG_TAIL = os.path.split(WRF_IMG)
                 if(WRF_IMG_TAIL.__str__().find(model)>=0 and WRF_IMG_TAIL.__str__().find(wrf_Data_Type)>=0):
                     Images.append(Image.open(WRF_SDV_Folder + "\\"+ WRF_IMG))
             for GCM_IMG in os.listdir(GCM_SDV_Folder): 
                 GCM_IMG_HEAD , GCM_IMG_TAIL = os.path.split(GCM_IMG)
                 if(GCM_IMG_TAIL.__str__().find(model)>=0 and GCM_IMG_TAIL.__str__().find(precipitation_Type)>=0):
                     Images.append(Image.open(GCM_SDV_Folder + "\\"+GCM_IMG))
             for Delta_IMG in os.listdir(Delta_SDV_Folder): 
                 Delta_IMG_HEAD , Delta_IMG_TAIL = os.path.split(Delta_IMG)
                 if(Delta_IMG_TAIL.__str__().find(model)>=0 and Delta_IMG_TAIL.__str__().find(precipitation_Type)>=0):
                     Images.append(Image.open(Delta_SDV_Folder + "\\"+Delta_IMG))
        
     
             widths, heights = zip(*(i.size for i in Images))
            
             total_width = sum(widths)
             max_height = max(heights)
            
             new_im = Image.new('RGBA', (total_width, max_height))
            
             x_offset = 0
             for im in Images:
                 new_im.paste(im, (x_offset,0))
                 x_offset += im.size[0]

             new_im.save("Ensemble SDV Plots"+ "\\" + "Ensemble SDV Plot with "+ model + " based on " + precipitation_Type+".png")


#stitch TOE and delta plots for each model pair together
def Generate_Three_Pane_TOE_Plots():
    models = ["access1.3","bcc","canesm2","ccsm4","csiro","fgoals","gfdl","giss","miroc5","mri","noresm1"]
    WRF_TOE_Folder = r"WRF TOE Plots" 
    GCM_TOE_Folder = r"GCM TOE Plots"
    Delta_TOE_Folder = r"ToE Deltas"
    
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
             for WRF_IMG in os.listdir(WRF_TOE_Folder):
                 WRF_IMG_HEAD , WRF_IMG_TAIL = os.path.split(WRF_IMG)
                 if(WRF_IMG_TAIL.__str__().find(model)>=0 and WRF_IMG_TAIL.__str__().find(wrf_Data_Type)>=0):
                     Images.append(Image.open(WRF_TOE_Folder + "\\"+ WRF_IMG))
             for GCM_IMG in os.listdir(GCM_TOE_Folder): 
                 GCM_IMG_HEAD , GCM_IMG_TAIL = os.path.split(GCM_IMG)
                 if(GCM_IMG_TAIL.__str__().find(model)>=0 and GCM_IMG_TAIL.__str__().find(temperature_Type)>=0):
                     Images.append(Image.open(GCM_TOE_Folder + "\\"+GCM_IMG))
             for Delta_IMG in os.listdir(Delta_TOE_Folder): 
                 Delta_IMG_HEAD , Delta_IMG_TAIL = os.path.split(Delta_IMG)
                 if(Delta_IMG_TAIL.__str__().find(model)>=0 and Delta_IMG_TAIL.__str__().find(temperature_Type)>=0):
                     Images.append(Image.open(Delta_TOE_Folder + "\\"+Delta_IMG))
        
     
             widths, heights = zip(*(i.size for i in Images))
            
             total_width = sum(widths)
             max_height = max(heights)
            
             new_im = Image.new('RGBA', (total_width, max_height))
            
             x_offset = 0
             for im in Images:
                 new_im.paste(im, (x_offset,0))
                 x_offset += im.size[0]

             new_im.save("Ensemble TOE Plots"+ "\\" + "Ensemble TOE Plot with " + model + " based on " + temperature_Type+".png")
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
             for WRF_IMG in os.listdir(WRF_TOE_Folder):
                 WRF_IMG_HEAD , WRF_IMG_TAIL = os.path.split(WRF_IMG)
                 if(WRF_IMG_TAIL.__str__().find(model)>=0 and WRF_IMG_TAIL.__str__().find(wrf_Data_Type)>=0):
                     Images.append(Image.open(WRF_TOE_Folder + "\\"+ WRF_IMG))
             for GCM_IMG in os.listdir(GCM_TOE_Folder): 
                 GCM_IMG_HEAD , GCM_IMG_TAIL = os.path.split(GCM_IMG)
                 if(GCM_IMG_TAIL.__str__().find(model)>=0 and GCM_IMG_TAIL.__str__().find(precipitation_Type)>=0):
                     Images.append(Image.open(GCM_TOE_Folder + "\\"+GCM_IMG))
             for Delta_IMG in os.listdir(Delta_TOE_Folder): 
                 Delta_IMG_HEAD , Delta_IMG_TAIL = os.path.split(Delta_IMG)
                 if(Delta_IMG_TAIL.__str__().find(model)>=0 and Delta_IMG_TAIL.__str__().find(precipitation_Type)>=0):
                     Images.append(Image.open(Delta_TOE_Folder + "\\"+Delta_IMG))
        
     
             widths, heights = zip(*(i.size for i in Images))
            
             total_width = sum(widths)
             max_height = max(heights)
            
             new_im = Image.new('RGBA', (total_width, max_height))
            
             x_offset = 0
             for im in Images:
                 new_im.paste(im, (x_offset,0))
                 x_offset += im.size[0]

             new_im.save("Ensemble TOE Plots"+ "\\" + "Ensemble TOE Plot with "+ model + " based on " + precipitation_Type+".png")
    
        
#execute
Generate_Three_Pane_TOE_Plots()
Generate_Three_Pane_SDV_Plots()
Generate_Three_Pane_Slope_Plots()