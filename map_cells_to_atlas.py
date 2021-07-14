#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 16:18:47 2021

@author: kepecs

map transformed cells to region
"""

import os, tifffile as tif, json, pandas as pd, numpy as np, seaborn as sns
import ClearMap.Alignment.Annotation as ano     
from collections import Counter
import matplotlib.pyplot as plt, matplotlib as mpl
from scipy.stats import ttest_ind

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6

def findkeys(node, kv):
    """
    iterate through allen hierarchy in json file
    """
    if isinstance(node, list):
        for i in node:
            for x in findkeys(i, kv):
               yield x
    elif isinstance(node, dict):
        if kv in node:
            yield node[kv]
        for j in node.values():
            for x in findkeys(j, kv):
                yield x

def get_progeny(dic,parent_structure,progeny_list):
    """austin hoag's function to get progeny from allen hierarchy"""
    if "msg" in list(dic.keys()): 
        dic = dic["msg"][0]
    
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

src = "/home/kepecs/Documents/"
annpth = "/home/kepecs/python/ClearMap2/ClearMap/Resources/Atlas/ABA_25um_annotation.tif"
ontology_file = "/home/kepecs/python/ClearMap2/ClearMap/Resources/Atlas/ABA_annotation.json"
animals = ["AA6-AK1a", "AA6-AK1b", "AA6-AK1c", "AA6-AK1d",
           "AA6-AK3b", "AA6-AK3c"]
#make annotation df
df = pd.DataFrame()
ann = tif.imread(annpth)
with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)
ontology_dict = ontology_dict["msg"][0]    
anndf = pd.DataFrame()
for k,v in ontology_dict.items():    
    if k != "children": #ignore children ids, can get from json
        anndf[k] = list(findkeys(ontology_dict, k))
#get total counts per region
anndf["total_voxels"] = [ann[ann==iid].shape[0] for iid in anndf.id.values]
df = anndf.copy()
df["cell_count"] = np.zeros(len(anndf))
animaldfs = []
#iterate thru the animals
for animal in animals:
    print(animal)
    pth = os.path.join(src, animal)
    animaldf = df.copy()
    animaldf["animal"] = animal
    #561
    pnts = os.path.join(pth, "{0}_561_points/posttransformed_zyx_voxels.npy".format(animal))
    pnts = np.load(pnts)
    labels = []
    for pnt in pnts:
        pnt = pnt.astype(int)
        try:
            labels.append(ann[pnt[0],pnt[1],pnt[2]])
        except Exception as e:
            print(e)
    counts = Counter(labels)
    for k,v in counts.items():
        animaldf.loc[animaldf.id == k, "cell_count"] = v #assign counts
    animaldfs.append(animaldf)
#concat
animaldfmaster561 = pd.concat(animaldfs)
#export 
animaldfmaster561.to_csv("/home/kepecs/Documents/561_cell_counts_all_animals.csv", index = False)   
#640
animaldfs = []
for animal in animals:
    print(animal)
    pth = os.path.join(src, animal)
    animaldf = df.copy()
    animaldf["animal"] = animal
    pnts = os.path.join(pth, "{0}_640_points/posttransformed_zyx_voxels.npy".format(animal))
    #reorient to xyz
    pnts = np.load(pnts) #get cells detected
    labels = []
    for pnt in pnts:
        pnt = pnt.astype(int)
        try:
            labels.append(ann[pnt[0],pnt[1],pnt[2]])
        except Exception as e:
            print(e)
    counts = Counter(labels)
    for k,v in counts.items():
        animaldf.loc[animaldf.id == k, "cell_count"] = v #assign counts
    animaldfs.append(animaldf)
#concat
animaldfmaster640 = pd.concat(animaldfs)
#export 
animaldfmaster640.to_csv("/home/kepecs/Documents/640_cell_counts_all_animals.csv", index = False)     

#%% 
#inspect
#assign groups
#scale
scale = 0.025 #mm
testan = ["AA6-AK1a", "AA6-AK1b", "AA6-AK1c", "AA6-AK1d"]
animaldfmaster561.loc[animaldfmaster561.animal.isin(testan), "group"] = "punished"
animaldfmaster561.loc[~animaldfmaster561.animal.isin(testan), "group"] = "control"
#filter out 0 vox regions
analyse561 = animaldfmaster561.copy()
analyse561 = analyse561[analyse561.total_voxels > 0]
analyse561["density"] = analyse561["cell_count"]/(analyse561["total_voxels"]*(scale**3))
#normalize by animal
for animal in animals:
    analyse561.loc[analyse561.animal == animal, "frac_of_cell_count"] = analyse561.loc[analyse561.animal == animal, 
                                "cell_count"]/analyse561.loc[analyse561.animal == animal, "cell_count"].sum()
#%cells/mm3
analyse561["frac_of_density"] = analyse561["frac_of_cell_count"]/(analyse561["total_voxels"]*(scale**3))

iid_pval = {}
#density p-values
for iid in np.unique(animaldfmaster561.id):
    tstat,pval = ttest_ind(analyse561.loc[(analyse561.id == iid) & 
                           (analyse561.group == "punished"), "frac_of_density"].values,
                             analyse561.loc[(analyse561.id == iid) & 
                           (analyse561.group == "control"), "frac_of_density"].values, equal_var = False)
    #assign
    analyse561.loc[(analyse561.id == iid) , "uncorrected p-value (frac_of_density)"] = pval
    if str(pval) != "nan": iid_pval[iid] = pval
#correct...
from statsmodels.stats.multitest import multipletests
pvals = np.array(list(iid_pval.values()))
r,pvals_corr,alphas,alphab = multipletests(pvals, method = "fdr_bh", alpha = 0.1)
for i,iid in enumerate(iid_pval.keys()):
    analyse561.loc[(analyse561.id == iid) , "fdr bh corrected p-value (frac_of_density)"] = pvals_corr[i]
#frac p-values
iid_pval = {}
for iid in np.unique(animaldfmaster561.id):
    tstat,pval = ttest_ind(analyse561.loc[(analyse561.id == iid) & 
                           (analyse561.group == "punished"), "frac_of_cell_count"].values,
                             analyse561.loc[(analyse561.id == iid) & 
                           (analyse561.group == "control"), "frac_of_cell_count"].values, equal_var = False)
    #assign
    analyse561.loc[(analyse561.id == iid) , "uncorrected p-value (frac_of_cell_count)"] = pval
    if str(pval) != "nan": iid_pval[iid] = pval
#correct...
from statsmodels.stats.multitest import multipletests
pvals = np.array(list(iid_pval.values()))
r,pvals_corr,alphas,alphab = multipletests(pvals, method = "fdr_bh", alpha = 0.1)
for i,iid in enumerate(iid_pval.keys()):
    analyse561.loc[(analyse561.id == iid) , "fdr bh corrected p-value (frac_of_cell_count)"] = pvals_corr[i]
analyse561.to_csv("/home/kepecs/Documents/561_cell_counts_all_animals_w_pvals.csv", index = False)   


#repeat for 640
#assign groups
#scale
scale = 0.025 #mm
testan = ["AA6-AK1a", "AA6-AK1b", "AA6-AK1c", "AA6-AK1d"]
animaldfmaster640.loc[animaldfmaster640.animal.isin(testan), "group"] = "punished"
animaldfmaster640.loc[~animaldfmaster640.animal.isin(testan), "group"] = "control"
#filter out 0 vox regions
analyse640 = animaldfmaster640.copy()
analyse640 = analyse640[animaldfmaster640.total_voxels > 0]
analyse640["density"] = analyse640["cell_count"]/(analyse640["total_voxels"]*(scale**3))
#normalize by animal
for animal in animals:
    analyse640.loc[analyse640.animal == animal, "frac_of_cell_count"] = analyse640.loc[analyse640.animal == animal, 
                                "cell_count"]/analyse640.loc[analyse640.animal == animal, "cell_count"].sum()
#%cells/mm3
analyse640["frac_of_density"] = analyse640["frac_of_cell_count"]/(analyse640["total_voxels"]*(scale**3))
   
#frac p-values
iid_pval = {}
for iid in np.unique(animaldfmaster561.id):
    tstat,pval = ttest_ind(analyse640.loc[(analyse640.id == iid) & 
                           (analyse640.group == "punished"), "frac_of_cell_count"].values,
                             analyse640.loc[(analyse640.id == iid) & 
                           (analyse640.group == "control"), "frac_of_cell_count"].values, equal_var = False)
    #assign
    analyse640.loc[(analyse640.id == iid) , "uncorrected p-value (frac_of_cell_count)"] = pval
    if str(pval) != "nan": iid_pval[iid] = pval
#correct...
from statsmodels.stats.multitest import multipletests
pvals = np.array(list(iid_pval.values()))
r,pvals_corr,alphas,alphab = multipletests(pvals, method = "fdr_bh", alpha = 0.1)
for i,iid in enumerate(iid_pval.keys()):
    analyse640.loc[(analyse640.id == iid) , "fdr bh corrected p-value (frac_of_cell_count)"] = pvals_corr[i]
    
#density frac p-values
iid_pval = {}
for iid in np.unique(animaldfmaster561.id):
    tstat,pval = ttest_ind(analyse640.loc[(analyse640.id == iid) & 
                           (analyse640.group == "punished"), "frac_of_density"].values,
                             analyse640.loc[(analyse640.id == iid) & 
                           (analyse640.group == "control"), "frac_of_density"].values, equal_var = False)
    #assign
    analyse640.loc[(analyse640.id == iid) , "uncorrected p-value (frac_of_density)"] = pval
    if str(pval) != "nan": iid_pval[iid] = pval
#correct...
from statsmodels.stats.multitest import multipletests
pvals = np.array(list(iid_pval.values()))
r,pvals_corr,alphas,alphab = multipletests(pvals, method = "fdr_bh", alpha = 0.1)
for i,iid in enumerate(iid_pval.keys()):
    analyse640.loc[(analyse640.id == iid) , "fdr bh corrected p-value (frac_of_density)"] = pvals_corr[i]
    
analyse640.to_csv("/home/kepecs/Documents/640_cell_counts_all_animals_w_pvals.csv", index = False)   

#%%
pval561 = analyse561[analyse561["uncorrected p-value (frac_of_cell_count)"]<0.05]
f, ax = plt.subplots(figsize=(7,14))
# g = sns.boxplot(x = "frac_of_cell_count", y = "name", hue = "group", data = pval561, orient = "h", showfliers = False, showcaps=False)
g = sns.stripplot(x = "frac_of_cell_count", y = "name", hue = "group", data = pval561, orient = "h", size = 4)
#hide the right and top spines
sns.despine(top=True, right=True, left=False, bottom=False)
g.set_xscale("symlog")
# g.set_xlim([-0.001, 0.02])
# ticks = np.arange(0,0.04,0.02)
# g.set_xticks(ticks); g.set_xticklabels(ticks)
ax.set_xlabel("cell count/total cells in brain (log)")
ax.set_ylabel("brain region")
ax.set_title("channel 561 nm, p < 0.05 (uncorrected)")
plt.savefig("/home/kepecs/Documents/561_boxplot_pval.svg", dpi = 300, bbox_inches="tight")

pval640 = analyse640[analyse640["uncorrected p-value (frac_of_cell_count)"]<0.05]
f, ax = plt.subplots(figsize=(7,17))
# g = sns.boxplot(x = "frac_of_cell_count", y = "name", hue = "group", data = pval640, orient = "h", showfliers = False, showcaps=False)
g = sns.stripplot(x = "frac_of_cell_count", y = "name", hue = "group", data = pval640, orient = "h", size = 4)
#hide the right and top spines
sns.despine(top=True, right=True, left=False, bottom=False)
g.set_xscale("symlog")
# g.set_xlim([-0.001, 0.04])
# ticks = np.arange(0,0.05,0.01)
# g.set_xticks(ticks); g.set_xticklabels(ticks)
g.set_yticklabels(g.get_yticklabels(), fontsize = 7)
ax.set_xlabel("cell count/total cells in brain (log)")
ax.set_ylabel("brain region")
ax.set_title("channel 640 nm, p < 0.05 (uncorrected)")

plt.savefig("/home/kepecs/Documents/640_boxplot_pval.svg", dpi = 300, bbox_inches="tight")

#%%
#get stripplot of all regions
sois = ["Isocortex", "Thalamus", "Hypothalamus", "Midbrain", "Hindbrain"]
ordered_sois = []
for soi in sois:
    get_progeny(ontology_dict, soi, ordered_sois)
#561    
#sort by nonzero sois
ordered_sois = [xx for xx in ordered_sois if any(analyse561.loc[analyse561.name == xx, "cell_count"] > 0)]
ordered_acronym = [analyse561.loc[analyse561.name == soi, "acronym"].values[0] for soi in ordered_sois]
    
f, ax = plt.subplots(figsize=(25,7))
g = sns.stripplot(x = "acronym", y = "frac_of_density", hue = "group", order = ordered_acronym,
                  data = analyse561, orient = "v", size = 4)
#shutoff ticks
g.set_xticklabels(ordered_acronym, rotation = 90, Fontsize=4)# g.set_ytickslabels()
# g.set_ylim([-1,5])
plt.savefig("/home/kepecs/Documents/561_manhattan_plot_frac_of_density.svg", dpi = 300, bbox_inches="tight")
#640
#sort by nonzero sois
ordered_sois = [xx for xx in ordered_sois if any(analyse640.loc[analyse561.name == xx, "cell_count"] > 0)]
ordered_acronym = [analyse640.loc[analyse640.name == soi, "acronym"].values[0] for soi in ordered_sois]
    
f, ax = plt.subplots(figsize=(25,7))
g = sns.stripplot(x = "acronym", y = "frac_of_density", hue = "group", order = ordered_acronym,
                  data = analyse640, orient = "v", size = 4)
#shutoff ticks
g.set_xticklabels(ordered_acronym, rotation = 90, Fontsize=4)# g.set_ytickslabels()
plt.savefig("/home/kepecs/Documents/640_manhattan_plot_frac_of_density.svg", dpi = 300, bbox_inches="tight")

