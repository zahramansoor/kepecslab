#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 16:40:44 2021

@author: kepecs
"""

import os, tifffile as tif, numpy as np

src = "/home/kepecs/Documents/density_vols"

test561 = [os.path.join(src, xx) for xx in os.listdir(src) if "AA6-AK1" in xx and "561" in xx]
control561 = [os.path.join(src, xx) for xx in os.listdir(src) if "AA6-AK3" in xx and "561" in xx]
test640 = [os.path.join(src, xx) for xx in os.listdir(src) if "AA6-AK1" in xx and "640" in xx]
control640 = [os.path.join(src, xx) for xx in os.listdir(src) if "AA6-AK3" in xx and "640" in xx]
#read imgs
imtest561 = np.array([tif.imread(xx) for xx in test561])
imcontrol561 = np.array([tif.imread(xx) for xx in control561])
imtest640 = np.array([tif.imread(xx) for xx in test640])
imcontrol640 = np.array([tif.imread(xx) for xx in control640])

imtest561_mean = np.mean(imtest561, axis=0)
imtest640_mean = np.mean(imtest640, axis=0)
imcontrol561_mean = np.mean(imcontrol561, axis=0)
imcontrol640_mean = np.mean(imcontrol640, axis=0)

tif.imsave("/home/kepecs/Documents/density_vols/test561mean.tif", imtest561_mean.astype("uint16"))
tif.imsave("/home/kepecs/Documents/density_vols/test640mean.tif", imtest640_mean.astype("uint16"))
tif.imsave("/home/kepecs/Documents/density_vols/control561mean.tif", imcontrol561_mean.astype("uint16"))
tif.imsave("/home/kepecs/Documents/density_vols/control640mean.tif", imcontrol640_mean.astype("uint16"))