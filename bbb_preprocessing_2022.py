#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 10:57:45 2022

@author: kepecs
"""

import tifffile as tif, matplotlib.pyplot as plt, numpy as np, os, pandas as pd, json, copy
from scipy.stats import ttest_ind as ttest
from statsmodels.stats.multitest import multipletests
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
#preprocessing to make into dataframe like amy's

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
    
#get voxel counts from brainpipe
df = pd.DataFrame()
allen = pd.read_excel("/home/kepecs/Documents/allen_id_table_w_voxel_counts.xlsx")
df["name"] = allen.name
df["voxels_in_structure"] = allen["voxels_in_structure"]

from allensdk.api.queries.ontologies_api import OntologiesApi
oapi = OntologiesApi()
structure_graph = oapi.get_structures_with_sets([1])  # 1 is the id of the adult mouse structure graph
from allensdk.core.structure_tree import StructureTree

# This removes some unused fields returned by the query
structure_graph = StructureTree.clean_structures(structure_graph)  
tree = StructureTree(structure_graph)
# get the ids of all the structure sets in the tree
structure_set_ids = tree.get_structure_sets()

# query the API for information on those structure sets
allen_stuff = pd.DataFrame(oapi.get_structure_sets(structure_set_ids))
target_group = 167587189 #summary structures
summary_structures = tree.get_structures_by_set_id([target_group])
#%%
sum_structs = pd.DataFrame(summary_structures)["name"]
#count voxels in regions
vox = pd.DataFrame()
vox["name"] = sum_structs
for area in sum_structs: 
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
    for soi in sum_structs:
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
#Calculate % counts
animals = ['cadw2_rhrh_ca', 'pr2w2_lhrh_ca', 'cadw1_lh_cach', 'pr2w2_rh_ctrl', 'pr2w2_rhrh_ct',
       'cadw2_nh_ctrl']

for nm in animals:
    total = np.nansum(df[nm])
    df[nm+"_percent_count"] = [xx/total if xx != np.nan else np.nan for xx in df[nm]]
#%%            
#filter out nan
import seaborn as sns
from matplotlib.colors import LogNorm

dfp = df.dropna(how = "all", subset = animals)
dfp.index = dfp.name
dfp = dfp.drop(columns = "name")
#add voxel column
for area in dfp.index:
    if vox.loc[vox.name == area, "voxels_in_structure"].values>0:
        dfp.loc[dfp.index == area, "voxels_in_structure"] = vox.loc[vox.name == area, "voxels_in_structure"].values

dfd = dfp.copy()
#get density
for animal in animals:
    dfd[animal] = dfp[animal]/dfp["voxels_in_structure"]
dfd = dfd.drop(columns = "voxels_in_structure").dropna(how = "all", subset = animals)[animals]   
#%% 
plt.figure(figsize=(30,1.5))
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
#shuffle regions
an = ['cadw2_rhrh_ca_percent_count', 'pr2w2_lhrh_ca_percent_count', 'cadw1_lh_cach_percent_count']
arr = np.array(dfp[an])
shufs = [arr[np.random.choice(np.arange(len(arr)), replace=False, size=len(arr)),:] for i in range(10000)]
shufmean = np.mean(shufs, axis=0)
#one sided ttest
# dfd["pvalue"]= [ttest(arr[i], shufmean[i])[1] if ttest(arr[i], shufmean[i])[0]>0 else 1 for i in range(len(arr))]
#two sided  ttest
dfp["pvalue"]= [ttest(arr[i], shufmean[i])[1] for i in range(len(arr))]
dfp["qvalue"] = multipletests(dfp.pvalue.values, method="fdr_bh")[1]
#only get regions which dont have all zeros
dfpp = dfp[(dfp[an[0]]!=0) | (dfp[an[1]]!=0) | (dfp[an[2]]!=0)]
dfpp.to_csv("/home/kepecs/Desktop/text.csv")

plt.figure(figsize=(30,1.5))
cmap = copy.copy(plt.cm.Blues)#plt.cm.Reds)
cmap.set_over(plt.cm.Blues(1.0)) #cmap.set_over('maroon')
cmap.set_under('w')

## WORK IN PROGRESS
p = sns.heatmap(dfpp.drop(columns = "voxels_in_structure").T, xticklabels = dfpp.drop(columns = "voxels_in_structure").index, cmap = cmap, 
                # norm = LogNorm(),
                # vmin=0, vmax=1e-2, 
                cbar_kws={'label': 'cadaverine + cells / total cadaverine + cells in brain'})
#how to quantify density?
p.set_xticklabels(dfpp.drop(columns = "voxels_in_structure").index, size = 8)
plt.savefig("/home/kepecs/Desktop/p_count.jpg", bbox_inches = "tight")

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
