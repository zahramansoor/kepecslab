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

#%%
from allensdk.core.reference_space_cache import ReferenceSpaceCache
from allensdk.api.queries.ontologies_api import OntologiesApi
oapi = OntologiesApi()
import pandas as pd
from allensdk.core.structure_tree import StructureTree

# This removes some unused fields returned by the query
structure_graph = oapi.get_structures_with_sets([1])  # 1 is the id of the adult mouse structure graph
structure_graph = StructureTree.clean_structures(structure_graph)  
tree = StructureTree(structure_graph)
# get the ids of all the structure sets in the tree
structure_set_ids = tree.get_structure_sets()

# query the API for information on those structure sets
allen_stuff = pd.DataFrame(oapi.get_structure_sets(structure_set_ids))
target_group = 167587189 #summary structures
summary_structures = tree.get_structures_by_set_id([target_group])

reference_space_key = 'annotation/ccf_2017'
resolution = 25
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
# ID 1 is the adult mouse structure graph
tree = rspc.get_structure_tree(structure_graph_id=1) 

structure_set_ids = tree.get_structure_sets()

# query the API for information on those structure sets
allen_stuff = pd.DataFrame(oapi.get_structure_sets(structure_set_ids))

import os
annotation, meta = rspc.get_annotation_volume()
# The file should appear in the reference space key directory
os.listdir(reference_space_key)
rsp = rspc.get_reference_space()

# This gets all of the structures targeted by the Allen Brain Observatory project
# brain_observatory_structures = rsp.structure_tree.get_structures_by_set_id([167587189])
# brain_observatory_ids = [st['name'] for st in brain_observatory_structures]
# brain_observatory_ids 

ventricular_mask = rsp.make_structure_mask([73])

# view in horizontal section
fig, ax = plt.subplots(figsize=(10, 10))
#look at ventricular mask
plt.imshow(np.max(ventricular_mask, axis=0), interpolation='none', cmap=plt.cm.afmhot)

#make mask of edge of the brain
brain_mask = rsp.make_structure_mask([8])
#dilate around mask
from skimage.morphology import binary_dilation
from skimage import morphology
from scipy.ndimage.morphology import distance_transform_edt
distance_space_outside = distance_transform_edt(np.logical_not(brain_mask.astype("bool")), 
                                                sampling=(25,25,25)) #INSIDE

mask = np.copy(distance_space_outside)
struct_microns_to_dilate = 80
mask[distance_space_outside >= struct_microns_to_dilate] = 0
outline = brain_mask-mask
#zero out other values <=1 not related to cdists
outline[outline>=1] = 0
#add to ventricular outlines
outline_w_ventricles = (outline.astype(bool).astype(int) + ventricular_mask).astype(bool).astype(int)
#%%
#get lists of structures used in statistics
structure_graph = oapi.get_structures_with_sets([1])  # 1 is the id of the adult mouse structure graph

structure_graph = StructureTree.clean_structures(structure_graph)  
tree = StructureTree(structure_graph)
# get the ids of all the structure sets in the tree
structure_set_ids = tree.get_structure_sets()

# query the API for information on those structure sets
allen_stuff = pd.DataFrame(oapi.get_structure_sets(structure_set_ids))
target_group = 167587189 #summary structures
summary_structures = tree.get_structures_by_set_id([target_group])

#get iids of structures I'm analysis
sum_structs = pd.DataFrame(summary_structures)["id"].values

import functools
from allensdk.core.reference_space import ReferenceSpace
# Define a wrapper function that will control the mask generation. 
# This one checks for a nrrd file in the specified base directory 
# and builds/writes the mask only if one does not exist
annotation_dir = '/home/kepecs/Documents/cadaverine_slices/annotation'
mask_writer = functools.partial(ReferenceSpace.check_and_write, annotation_dir)
    
rsp = rspc.get_reference_space()
# many_structure_masks is a generator - nothing has actrually been run yet
mask_generator = rsp.many_structure_masks(sum_structs, mask_writer)

# consume the resulting iterator to make and write the masks
for structure_id in mask_generator:
    print( 'made mask for structure {0}.'.format(structure_id) ) 

os.listdir(annotation_dir)
#%%
#calculate distance from ventricles
#calculate centroids
from scipy.ndimage.measurements import center_of_mass
from scipy.spatial.distance import cdist
import nrrd
annotation_dir = '/home/kepecs/Documents/cadaverine_slices/annotation'

#centroid df
sum_structs = pd.DataFrame(summary_structures)["id"].values
df = pd.DataFrame(); df["id"] = sum_structs; df["name"] = pd.DataFrame(summary_structures)["name"].values
#find nonzero points of ventricular mask
v = np.array([np.nonzero(outline_w_ventricles)[0],np.nonzero(outline_w_ventricles)[1],np.nonzero(outline_w_ventricles)[2]]).T
for mask in os.listdir(annotation_dir):
    readdata, header = nrrd.read(os.path.join(annotation_dir, mask))
    z, y, x = center_of_mass(readdata)
    df.loc[df.id == int(mask[10:-5]), "centroid_z"] = int(z)
    df.loc[df.id == int(mask[10:-5]), "centroid_y"] = int(y)
    df.loc[df.id == int(mask[10:-5]), "centroid_x"] = int(x)
    print(mask[10:-5])
    dists = np.array([np.linalg.norm(np.array([z,y,x])-v[i]) for i in range(len(v))])
    df.loc[df.id == int(mask[10:-5]), "median_euclidean_dist_from_ventricle_borders"] = np.median(dists)
    df.loc[df.id == int(mask[10:-5]), "min_euclidean_dist_from_ventricle_borders"] = np.min(dists)
    df.loc[df.id == int(mask[10:-5]), "max_euclidean_dist_from_ventricle_borders"] = np.max(dists)
#export    
df.to_csv("/home/kepecs/Documents/cadaverine_slices/distance_from_ventricle.csv", index = None)

    
