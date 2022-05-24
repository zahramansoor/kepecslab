#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 10:57:45 2022

@author: kepecs
"""

import tifffile as tif, matplotlib.pyplot as plt, numpy as np, os, pandas as pd, json, copy

#%%
src = "/home/kepecs/Documents/cadaverine_slices/original"
ims = [os.path.join(src, xx) for xx in os.listdir(src) if "done" not in xx and "ch00" in xx]
for im in ims[22:]:
    print(im)
    dst = "/home/kepecs/Documents/cadaverine_slices/modified_ch01"
    img = tif.imread(im)
    img = np.rot90(img, axes=(1,0))
    pad = 5000
    imgb = np.zeros((img.shape[0]+pad, img.shape[1]+pad))
    imgb[int(pad/2):imgb.shape[0]-int(pad/2), int(pad/2):imgb.shape[1]-int(pad/2)]=img
    plt.imshow(imgb)
    tif.imsave(os.path.join(dst, os.path.basename(im)), imgb.astype('uint16'))
    del(img); del(imgb)

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
    csv_soi = csv_["name"].value_counts() # cell counts
    #csv_soi = csv_["name"].value_counts()/len(csv_) #% cell counts
    for soi in csv_soi.index:
        if soi != "Basic cell groups and regions" and soi != "Cerebral cortex":
            try: #if an entry already exists
                df.loc[df.name == soi, nm] = df.loc[df.name == soi, nm].sum() + csv_soi[soi]
            except:
                df.loc[df.name == soi, nm] = csv_soi[soi]
#%%            
#filter out nan
import seaborn as sns
animals = ['cadw2_rhrh_cachexia', 'pr2w2_lhrh_cachexia_series', 'cadw2_lh_cachexia', 
           'pr2w2_rh_ctrl_serie', 'pr2w2_rhrh_ctrl_', ]
dfp = df.dropna(how = "all", subset = animals)
dfp.index = dfp.name
dfp = dfp.drop(columns = "name")
dfp = dfp[animals]

plt.figure(figsize=(1,18))
cmap = copy.copy(plt.cm.Blues)#plt.cm.Reds)
cmap.set_over(plt.cm.Blues(1.0)) #cmap.set_over('maroon')
cmap.set_under('w')
p = sns.heatmap(dfp, xticklabels = animals, yticklabels = dfp.index, cmap = cmap, vmin=0, vmax=750,
            cbar_kws={'label': 'cadaverine + cells'})
p.set_yticklabels(dfp.index, size = 7)
plt.savefig("/home/kepecs/Desktop/test.pdf", bbox_inches = "tight")