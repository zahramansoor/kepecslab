#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 10:44:47 2021

@author: kepecs
"""

import tifffile as tif, numpy as np, os, matplotlib.pyplot as plt, sys
import shutil, cv2 
from skimage.morphology import ball

src = "/home/kepecs/Documents/"
dst = "/home/kepecs/Documents/checks"
if not os.path.exists(dst): os.mkdir(dst)
atl_pth = "/home/kepecs/python/ClearMap2/ClearMap/Resources/Atlas/ABA_25um_reference.tif"
atl = tif.imread(atl_pth)
atlz,atly,atlx = atl.shape
animals = ["AA6-AK3b", "AA6-AK3c"]
# ["AA6-AK1d", #["AA6-AK1a", "AA6-AK1b", "AA6-AK1c", 
           # ]
chs = [561, 640]
for animal in animals:
    print(animal)
    for ch in chs:
        print(ch)
        converted_points = os.path.join(src, animal, "{0}_{1}_points".format(animal, ch), "posttransformed_zyx_voxels.npy")
        converted_points = np.load(converted_points)
        arr=converted_points.astype(int)
        cell=np.zeros((atlz,atly,atlx)) #init cellmap
        miss = 0
        for pnt in arr:
            z,y,x=pnt
            try:
                cell[z,y,x] = 1
            except:
                miss+=1
        #apply x y dilation
        r = 3
        selem = ball(r)[int(r/2)]
        cell_map = cell.astype("uint8")
        cell_map = np.asarray([cv2.dilate(cell_map[i], selem, iterations = 1) for i in range(cell_map.shape[0])])
        
        #overlay of atlas
        merged = np.stack([atl, cell_map, np.zeros_like(cell_map)], -1)
        tif.imsave(os.path.join(dst, "{0}_{1}_overlay.tif".format(animal, ch)), merged.astype("uint16"))

#%%
#animal name 
animal = "AA6-AK1d"
channel = 640
#raw point overlay
raw = "/home/kepecs/Documents/{0}/{0}_cells_{1}_filtered20.npy".format(animal, channel)
src = "/mnt/uncertainty/{0}/{0}_{1}".format(animal, channel)
zrng = [250,260]
pnts = np.load(raw)
pnts = pnts[(pnts["z"]>=zrng[0]) & (pnts["z"]<zrng[1])]
try:
    imgs = [os.path.join(src,xx) for xx in os.listdir(src) if "tif" in xx and int(xx[28:31])>=zrng[0] and int(xx[28:31])<zrng[1]]; imgs.sort()
except:
    imgs = [os.path.join(src,xx) for xx in os.listdir(src) if "tif" in xx and int(xx[27:30])>=zrng[0] and int(xx[27:30])<zrng[1]]; imgs.sort()
imgs = np.array([tif.imread(xx) for xx in imgs])
z,y,x = imgs.shape
cell=np.zeros((zrng[1]-zrng[0],y,x)) #init cellmap

for pnt in pnts:
    cell[pnt['z']-zrng[0],pnt['y'],pnt['x']] = 1
#apply x y dilation
r = 6
selem = ball(r)[int(r/2)]
cell_map = cell.astype("uint8")
cell_map = np.asarray([cv2.dilate(cell_map[i], selem, iterations = 1) for i in range(cell_map.shape[0])])
#overlay of cell map
merged = np.stack([np.max(imgs,axis=0), np.max(cell_map,axis=0), np.zeros((y,x))], -1)
tif.imsave("/home/kepecs/Documents/{0}_raw_{1}_z{2}-{3}.tif".format(animal, 
                channel, zrng[0],zrng[1]), merged.astype("uint16"))

#%%
import pandas as pd 

json = "/home/kepecs/python/ClearMap2/ClearMap/Resources/Atlas/ABA_annotation.json"
allen = pd.read_json(json)