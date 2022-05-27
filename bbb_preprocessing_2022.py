#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 10:57:45 2022

@author: kepecs
"""

import tifffile as tif, matplotlib.pyplot as plt, numpy as np, os, pandas as pd, json, copy

#%%
src = "/home/kepecs/Documents/cadaverine_slices/original"
ims = [os.path.join(src, xx) for xx in os.listdir(src) if "done" not in xx and "ch01" in xx and "cadw2_lh" in xx]
for im in ims:
    print(im)
    dst = "/home/kepecs/Documents/cadaverine_slices/modified_ch01"
    img = tif.imread(im)
    img = np.rot90(img, axes=(1,0))
    pad = 5000
    imgb = np.zeros((img.shape[0]+pad, img.shape[1]+pad))
    imgb[int(pad/2):imgb.shape[0]-int(pad/2), int(pad/2):imgb.shape[1]-int(pad/2)]=img
    tif.imsave(os.path.join(dst, os.path.basename(im)), imgb.astype('uint16'))
    del(img); del(imgb)

#%%
#postprocessing
fl = "/home/kepecs/Documents/cadaverine_slices/modified_ch01/csv_data"
ontology_file = "/home/kepecs/Documents/allen.json"
csvs = [os.path.join(fl, xx) for xx in os.listdir(fl) if "csv" in xx]

df = pd.DataFrame()
allen = pd.read_excel("/home/kepecs/Documents/allen_id_table_w_voxel_counts.xlsx")
df["name"] = allen.name
df["voxels_in_structure"] = allen["voxels_in_structure"]

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

#make a list of sois you want to quantify
areas = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Frontal pole, cerebral cortex", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area", "Entorhinal area", "Ventral posteromedial nucleus of the thalamus", "Ventral posterolateral nucleus of the thalamus",
            "Ventral anterior-lateral complex of the thalamus", "Anteroventral nucleus of thalamus", 
            "Lateral dorsal nucleus of thalamus", "Paraventricular nucleus of the thalamus", "Medial habenula",
            "Lateral posterior nucleus of the thalamus", "Posterior triangular thalamic nucleus", "Mediodorsal nucleus of thalamus",
            "Posterior complex of the thalamus","Ventral medial nucleus of the thalamus","Reticular nucleus of the thalamus",
            "Lateral hypothalamic area","Periventricular region","Zona incerta","Mammillary body","Anterior hypothalamic nucleus",
            "Periventricular zone","Lateral preoptic area","Posterior hypothalamic nucleus", "Medial preoptic area","Tuberal nucleus",
            "Ventromedial hypothalamic nucleus","Medial preoptic nucleus","Medial mammillary nucleus",
            "Dorsomedial nucleus of the hypothalamus","Arcuate hypothalamic nucleus","Supramammillary nucleus","Paraventricular hypothalamic nucleus",
            "Fields of Forel","Anteroventral periventricular nucleus","Periventricular hypothalamic nucleus, intermediate part",
            "Caudoputamen", "Nucleus accumbens","Olfactory tubercle", "Lateral septal nucleus", "Pallidum"
            "Medial amygdalar nucleus", "Central amygdalar nucleus","Anterior amygdalar area", "Septofimbrial nucleus",
            "Fundus of striatum", "Intercalated amygdalar nucleus","Septohippocampal nucleus", "Claustrum", "Endopiriform nucleus", "Lateral amygdalar nucleus",
            "Basolateral amygdalar nucleus", "Basomedial amygdalar nucleus", "Posterior amygdalar nucleus"]

#count voxels in regions
vox = pd.DataFrame()
vox["name"] = areas
for area in areas: 
    progeny = []; get_progeny(ontology_dict, area, progeny)
    try: 
        voxs = [allen.loc[allen.name == area, "voxels_in_structure"].values]
    except:
        voxs = []
    for prog in progeny:
        try:
            voxs.append(allen.loc[allen.name == prog, "voxels_in_structure"].values)
        except:
            voxs.append(0)
    vox.loc[vox.name == area, "voxels_in_structure"] = np.sum(voxs)
        
#compile csvs
for csv in csvs:
    nm = os.path.basename(csv)[:13]
    csv_ = pd.read_csv(csv)
    csv_soi = csv_["name"].value_counts() # cell counts
    # csv_soi = csv_["name"].value_counts()/len(csv_) #% cell counts
    for soi in areas:
        if soi != "Basic cell groups and regions" and soi != "Cerebral cortex":
            progeny = []; get_progeny(ontology_dict, soi, progeny)
            #try except statements are to make sure all sois and children are included
            try:
                counts = [csv_soi[soi]]             
            except:
                counts = []
            for prog in progeny:
                try: 
                    counts.append(csv_soi[prog])
                except:
                    counts.append(0)
            try: #if an entry already exists
                df.loc[df.name == soi, nm] = df.loc[df.name == soi, nm].sum() + np.sum(counts)
            except:
                df.loc[df.name == soi, nm] = np.sum(counts)
#%%            
#filter out nan
import seaborn as sns
from matplotlib.colors import LogNorm

animals = ['cadw2_rhrh_ca', 'pr2w2_lhrh_ca', 'cadw1_lh_cach', 'pr2w2_rh_ctrl', 'pr2w2_rhrh_ct',
       'cadw2_nh_ctrl']
dfp = df.dropna(how = "all", subset = animals)
dfp.index = dfp.name
dfp = dfp.drop(columns = "name")
dfp = dfp[animals]
#add voxel column
for area in areas:
    dfp.loc[dfp.index == area, "voxels_in_structure"] = vox.loc[vox.name == area, "voxels_in_structure"].values

dfd = dfp.copy()
#get density
for animal in animals:
    dfd[animal] = dfp[animal]/dfp["voxels_in_structure"]
dfd = dfd.drop(columns = "voxels_in_structure")    
plt.figure(figsize=(18,1.5))
cmap = copy.copy(plt.cm.Blues)#plt.cm.Reds)
cmap.set_over(plt.cm.Blues(1.0)) #cmap.set_over('maroon')
cmap.set_under('w')
p = sns.heatmap(dfd.T, yticklabels = animals, xticklabels = dfd.index, cmap = cmap, 
                # norm = LogNorm(),
                vmin=0, vmax=1e-2, 
                cbar_kws={'label': 'cadaverine + cells / total voxels in structure'})
#how to quantify density?
p.set_xticklabels(dfd.index, size = 8)
plt.savefig("/home/kepecs/Desktop/density.jpg", bbox_inches = "tight")
#%%
#split plots by regions
nc = ["Caudoputamen", "Nucleus accumbens","Olfactory tubercle", "Lateral septal nucleus", "Pallidum"]
plt.figure(figsize=(2,5))
cmap = copy.copy(plt.cm.Blues)#plt.cm.Reds)
cmap.set_over(plt.cm.Blues(1.0)) #cmap.set_over('maroon')
cmap.set_under('w')
p = sns.heatmap(dfp[dfd.index.isin(nc)].drop(columns = "voxels_in_structure"), xticklabels = animals, cmap = cmap, 
                # norm = LogNorm(),
                vmin=0, vmax=200, 
                cbar_kws={'label': 'cadaverine + cells'})
#how to quantify density?
p.set_yticklabels(p.get_yticklabels(), size = 8)

plt.savefig("/home/kepecs/Desktop/str_counts.jpg", bbox_inches = "tight")
