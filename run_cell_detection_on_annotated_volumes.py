#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 11:23:16 2021

@author: kepecs
"""

#ClearMap path
from scipy.io import loadmat
import sys, tifffile as tif, os, itertools, shutil, numpy as np
import argparse
sys.path.append('/home/kepecs/python/ClearMap2/')
from ClearMap.Environment import *  #analysis:ignore

volsrc = "/home/kepecs/Documents/cfos_annotations/volumes"
vols = [os.path.join(volsrc, xx) for xx in os.listdir(volsrc) if "tif" in xx]; vols.sort()

for vol in vols:
    volnm = os.path.basename(vol)[:-4]
    #make dest dir
    dst = os.path.join(volsrc,volnm)
    if not os.path.exists(dst): os.mkdir(dst)
    ch561 = [True if "561" in volnm else False][0]
    #parameter sweep
    bkshp = [None, (3,3), (5,5), (7,7), (10,10),
                   (13,13)]
    if ch561:
        shpthres = [1900,2100,2300,2600,3000] #diff params for 561 based on first param sweep
    else: shpthres = [100, 300, 500, 700, 900]
    maximashp = [5,10,15,20]
    # calculate number of iterations
    tick = 0
    for b, s, m in itertools.product(bkshp, shpthres, maximashp):
        tick +=1
    sys.stdout.write("\n\nNumber of iterations is {}:".format(tick))
    
    for i in range(tick):
        bkshp_, shpthres_, maximashp_ = [xx for xx in itertools.product(bkshp, shpthres, maximashp)][i] #parse out combinations
        print("\n")
        print("   iteration: {0}\n   background corr: {1}\n   shape detection: {2}".format(i,bkshp, shpthres_))
        print("\n")
        cell_detection_parameter = cells.default_cell_detection_parameter.copy();
        cell_detection_parameter["illumination_correction"] = None;
        if bkshp_ is not None: 
            cell_detection_parameter["background_correction"] = {"shape": bkshp_, "form": "Disk"}; 
        else: cell_detection_parameter["background_correction"] = bkshp_
        cell_detection_parameter["intensity_detection"]["measure"] = ["source"];
        cell_detection_parameter["shape_detection"]["threshold"] = shpthres_ #for 640 ch 500
        cell_detection_parameter["maxima_detection"]["shape"] = maximashp_
        
        processing_parameter = cells.default_cell_detection_processing_parameter.copy();
        processing_parameter.update(
               processes = 1, # 'serial', #multiple processes don't work on kepecs desktop bc of memory
               size_max = 20, #100, #35,
               size_min = 10,# 30, #30,
               optimization = False,
               optimization_fix = None,
               overlap  = 5, #32, #10,
               verbose = True
               )
        dstch = os.path.join(dst, "cells_bkshp{0}_shpthres{1}_maxshp{2}_raw.npy".format(bkshp_, shpthres_, maximashp_)) 
        cells.detect_cells(vol, dstch,
                             cell_detection_parameter=cell_detection_parameter, 
                             processing_parameter=processing_parameter)