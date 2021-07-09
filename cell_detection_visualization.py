#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 10:44:47 2021

@author: kepecs
"""

from scipy.io import loadmat
import tifffile as tif, numpy as np, os, matplotlib.pyplot as plt, sys
import shutil, cv2 
from skimage.morphology import ball

atl_pth = "/home/kepecs/python/ClearMap2/ClearMap/Resources/Atlas/ABA_25um_reference.tif"
atl = tif.imread(atl_pth)
z,y,x = atl.shape
#check
if isinstance(converted_points, str):
    converted_points = np.load(converted_points)
arr=converted_points.astype(int)
cell=np.zeros((z,y,x)) #init cellmap
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
plt.imshow(cell[300])
tif.imsave("/home/kepecs/Documents/test23.tif", merged)

#%%
#raw point overlay
raw = "/home/kepecs/Documents/AA6-AK1a/AA6-AK1a_cells_640_filtered20.npy"
src = "/mnt/uncertainty/AA6-AK1a/AA6-AK1a_640"
zrng = [600,630]
pnts = np.load(raw)
pnts = pnts[(pnts["z"]>=zrng[0]) & (pnts["z"]<zrng[1])]
imgs = [os.path.join(src,xx) for xx in os.listdir(src) if int(xx[28:31])>=zrng[0] and int(xx[28:31])<zrng[1]]; imgs.sort()
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
tif.imsave("/home/kepecs/Documents/raw_640_z{0}-{1}.tif".format(zrng[0],zrng[1]), merged.astype("uint16"))

#%%
import pandas as pd

json = "/home/kepecs/python/ClearMap2/ClearMap/Resources/Atlas/ABA_annotation.json"
allen = pd.read_json(json)