#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 17:07:08 2021

@author: kepecs
"""

import os, tifffile as tif, matplotlib.pyplot as plt, numpy as np, pandas as pd, seaborn as sns

fld = "/home/kepecs/Documents/slices/tifs/iba1/"

# imgs = [os.path.join(fld, xx) for xx in os.listdir(fld)]
# im = tif.imread(imgs[0])
# #display 
# plt.imshow(im[0, 6000:8000, 6000:8000]*70, cmap = "gist_yarg")
# #export
# for impth in imgs:
#     print(impth)
#     im = tif.imread(impth)
#     #channels
#     tif.imsave(os.path.join(fld, impth[:-4]+"_ch0.tif"), im[0].astype("uint16"))
#     tif.imsave(os.path.join(fld, impth[:-4]+"_ch1.tif"), im[1].astype("uint16"))
#     tif.imsave(os.path.join(fld, impth[:-4]+"_ch2.tif"), im[2].astype("uint16"))
#     tif.imsave(os.path.join(fld, impth[:-4]+"_ch3.tif"), im[3].astype("uint16"))
    
#%%
#get mean fluorescence across each slice for glia
#ch0 = dapi, ch1 = claudin, ch2 = glia, ch3 = cd31
ctrl_slices = [os.path.join(fld, xx) for xx in os.listdir(fld) if "ch2" in xx and "w2nh_nac" in xx]
test_slices = [os.path.join(fld, xx) for xx in os.listdir(fld) if "ch2" in xx and "w1lhrh_nac" in xx]
ctrl_slices = [tif.imread(xx) for xx in ctrl_slices]
test_slices = [tif.imread(xx) for xx in test_slices]
df = pd.DataFrame()
df["control_mean"] = [np.mean(xx) for xx in ctrl_slices] 
df["test_mean"] = [np.mean(xx) for xx in test_slices] 
df["control_integrated"] = [np.sum(xx) for xx in ctrl_slices] 
df["test_integrated"] = [np.sum(xx) for xx in test_slices] 

#%%

w1lhrh_nac = ["/home/kepecs/Documents/slices/10xtile_iba1_cd31_pilot_w1lhrh_nac_slice2_10xtile_corr_merged_ch2.csv",
              "/home/kepecs/Documents/slices/10xtile_iba1_cd31_pilot_w1lhrh_nac_slice3_10xtile_corr_merged_ch2.csv",
              "/home/kepecs/Documents/slices/10xtile_iba1_cd31_pilot_w1lhrh_nac_slice4_10xtile_corr_merged_ch2.csv"]
w2nh_nac = "/home/kepecs/Documents/slices/10xtile_iba1_cd31_pilot_w2nh_nac_slice2_10xtile_corr_merged_ch2.csv"

w1lhrh_nac = [pd.read_csv(xx) for xx in w1lhrh_nac]
w1lhrh_nac = pd.concat(w1lhrh_nac)
w2nh_nac = pd.read_csv(w2nh_nac)
w1lhrh_nac["group"] = ["cachexia"]*len(w1lhrh_nac)
w2nh_nac["group"] = ["control"]*len(w2nh_nac)
df = pd.concat([w2nh_nac,w1lhrh_nac])
sns.boxplot(x = "group", y = "intensity", data = df[df.name == "fiber tracts"], showfliers = False)
plt.figure()
sns.boxplot(x = "group", y = "intensity", data = df[df.name == "Nucleus accumbens"],showfliers = False)

testdf = df[df.name == "fiber tracts"]
from scipy.stats import ttest_ind
ttest_ind(testdf.loc[testdf.group == "control", "intensity"].values, testdf.loc[testdf.group == "cachexia", "intensity"].values)
