#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 10:57:45 2022

@author: kepecs
"""

import tifffile as tif, matplotlib.pyplot as plt, numpy as np, os, pandas as pd, json, copy

#%%
src = "/home/kepecs/Documents/cadaverine_slices"
ims = [os.path.join(src, xx) for xx in os.listdir(src) if "ch01" in xx and "pr2w2" in xx]
for im in ims:
    print(im)
    dst = "/home/kepecs/Documents/cadaverine_slices/modified_ch01"
    img = tif.imread(im)
    img = np.rot90(img, axes=(1,0))
    pad = 5000
    imgb = np.zeros((img.shape[0]+pad, img.shape[1]+pad))
    imgb[int(pad/2):imgb.shape[0]-int(pad/2), int(pad/2):imgb.shape[1]-int(pad/2)]=img
    plt.imshow(imgb)
    tif.imsave(os.path.join(dst, os.path.basename(im)), imgb.astype('uint16'))

#%%
#postprocessing
fl = "/home/kepecs/Documents/cadaverine_slices/modified_ch01"
ontology_file = "/home/kepecs/Documents/allen.json"
csvs = [os.path.join(fl, xx) for xx in os.listdir(fl) if "csv" in xx]

df = pd.DataFrame()
allen = pd.read_excel("/home/kepecs/Documents/allen_id_table_w_voxel_counts.xlsx")
df["name"] = allen.name

def get_progeny(dic,parent_structure,progeny_list):
   
    if "msg" in list(dic.keys()): dic = dic["msg"][0]
    
    name = dic.get("name")
    children = dic.get("children")
    if name == parent_structure:
        for child in children: # child is a dict
            child_name = child.get("name")
            progeny_list.append(child_name)
            get_progeny(child,parent_structure=child_name,progeny_list=progeny_list)
        return
    
    for child in children:
        child_name = child.get("name")
        get_progeny(child,parent_structure=parent_structure,progeny_list=progeny_list)
    return 

with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)
    
#compile csvs
for csv in csvs:
    nm = os.path.basename(csv)[:-18]
    csv_ = pd.read_csv(csv)
    csv_soi = csv_["name"].value_counts()/len(csv_) #% cell counts
    for soi in csv_soi.index:
        try: #if an entry already exists
            df.loc[df.name == soi, nm] = df.loc[df.name == soi, nm].sum() + csv_soi[soi]
        except:
            df.loc[df.name == soi, nm] = csv_soi[soi]
#%%            
#filter out nan
import seaborn as sns
dfp = df.dropna(subset = ["cadw2_rhrh_cachexia", "cadw2_lh_cachexia", "pr2w2_lhrh_cachexia_series"])
dfp.index = dfp.name
dfp = dfp.drop(columns = "name")

plt.figure(figsize=(1,18))
cmap = copy.copy(plt.cm.Blues)#plt.cm.Reds)
cmap.set_over(plt.cm.Blues(1.0)) #cmap.set_over('maroon')
cmap.set_under('w')
p = sns.heatmap(dfp, yticklabels = dfp.index, cmap = cmap, vmin=0, vmax=0.5,
            cbar_kws={'label': '% count'})
p.set_yticklabels(dfp.index, size = 8)