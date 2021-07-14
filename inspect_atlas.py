#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 09:49:08 2021

@author: kepecs
"""

import tifffile as tif, numpy as np, matplotlib.pyplot as plt
from collections import Counter

atlpth = "/home/kepecs/python/ClearMap2/ClearMap/Resources/Atlas/ABA_25um_annotation.tif"
atl = tif.imread(atlpth)
#parabrachial nuc ids
iids = [860, 868, 890]
#find indices correspond to the id
z,y,x=np.where(atl==iids[0])

plt.imshow(atl[z,y,x])
