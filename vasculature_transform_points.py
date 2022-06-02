#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 14:21:26 2022

@author: kepecs
Extract vasculature coordinates from Erturk's segmentation and export w/ registeration

"""
from scipy.io import loadmat
import tifffile as tif, numpy as np, os, matplotlib.pyplot as plt, sys
import shutil, cv2 
from skimage.morphology import ball
import pandas as pd
import json

#get cells
pth = "/home/kepecs/Documents/cadaverine_slices/downsized_for_atlas.tif"
img = tif.imread(pth)
z,y,x = np.nonzero(img)
arr = np.array([z,y,x]).T

def transformix_command_line_call(src, dst, transformfile):
    '''Wrapper Function to call transformix using the commandline, this can be time consuming
    
    Inputs
    -------------------
    src = volume path for transformation
    dst = folder to save file
    transformfile = final transform file from elastix registration
    
    '''
    from subprocess import check_output
    print ('Running transformix, this can take some time....\n')
    #sp.call(['transformix', '-in', src, '-out', dst, '-tp', transformfile])
    call = 'transformix -in {} -out {} -tp {}'.format(src, dst, transformfile)
    print(check_output(call, shell=True))
    print('Past transformix command line Call')      
            
    return

def change_transform_parameter_initial_transform(fl, initialtrnsformpth):
    '''
    (InitialTransformParametersFileName "NoInitialTransform")
    initialtrnsformpth = 'NoInitialTransform' or 'pth/to/transform.0.txt'
    '''
    fl1 = fl[:-5]+'0000.txt'
    with open(fl, 'r') as f, open(fl1, 'w') as new:
            for line in f:
                new.write('(InitialTransformParametersFileName "{}")\n'.format(initialtrnsformpth)) if 'InitialTransformParametersFileName' in line else new.write(line)
    #overwrite original transform file
    shutil.move(fl1, fl)
    return

def modify_transform_files(transformfiles, dst):
    """Function to copy over transform files, modify paths in case registration was done on the cluster, and tether them together
    
        Inputs
    ---------
    transformfiles = 
        list of all elastix transform files used, and in order of the original transform****
    """
    #new
    ntransformfiles = [os.path.join(dst, "order{}_{}".format(i,os.path.basename(xx))) for i,xx in enumerate(transformfiles)]    
    #copy files over
    [shutil.copy(xx, ntransformfiles[i]) for i,xx in enumerate(transformfiles)]   
    #modify each with the path
    for i,pth in enumerate(ntransformfiles):
        #skip first
        if i!=0:
            #read
            with open(pth, "r") as fl:
                lines = fl.readlines()
                fl.close()
            #copy
            nlines = lines
            #iterate
            for ii, line in enumerate(lines):
                if "(InitialTransformParametersFileName" in line:
                    nlines[ii] = "(InitialTransformParametersFileName {})\n".format(ntransformfiles[i-1])
            #write
            with open(pth, "w") as fl:
                for nline in lines:
                    fl.write(str(nline))
                fl.close()
    return ntransformfiles

def point_transformix(pretransform_text_file, transformfile, dst):
    """apply elastix transform to points      
    Inputs
    -------------
    pretransform_text_file = list of points that already have resizing transform
    transformfile = elastix transform file
    dst = folder
    
    Returns
    ---------------
    trnsfrm_out_file = pth to file containing post transformix points
    
    """
    sys.stdout.write("\n***********Starting Transformix***********")
    from subprocess import check_output
    #set paths    
    trnsfrm_out_file = os.path.join(dst, "outputpoints.txt")
    #run transformix point transform
    call = "transformix -def {} -out {} -tp {}".format(pretransform_text_file, dst, transformfile)
    print(check_output(call, shell=True))
    sys.stdout.write("\n   Transformix File Generated: {}".format(trnsfrm_out_file)); sys.stdout.flush()
    return trnsfrm_out_file

def create_text_file_for_elastix(src, dst):
    """
    Inputs
    ---------
    src = numpy file consiting of nx3 (ZYX points)
    dst = folder location to write points
    """
    print("This function assumes ZYX centers...")
    #setup folder
    if not os.path.exists(dst): os.mkdir(dst)
    #create txt file, with elastix header, then populate points
    pth=os.path.join(dst, "zyx_points_pretransform.txt")
    #load
    if type(src) == np.ndarray:
        arr = src
    else:
        arr = np.load(src) if src[-3:] == "npy" else loadmat(src)["cell_centers_orig_coord"]
    #convert
    stringtowrite = "\n".join(["\n".join(["{} {} {}".format(i[2], i[1], i[0])]) for i in arr]) ####this step converts from zyx to xyz*****
    #write file
    sys.stdout.write("writing centers to transfomix input points text file..."); sys.stdout.flush()
    with open(pth, "w+") as fl:
        fl.write("index\n{}\n".format(len(arr)))    
        fl.write(stringtowrite)
        fl.close()
    sys.stdout.write("...done writing centers\n"); sys.stdout.flush()
    return pth

def unpack_pnts(points_file, dst):
    """
    function to take elastix point transform file and return anatomical locations of those points
    
    Here elastix uses the xyz convention rather than the zyx numpy convention
    
    Inputs
    -----------
    points_file = post_transformed file, XYZ
    
    Returns
    -----------
    dst_fl = path to numpy array, ZYX
    
    """   
    #####inputs 
    assert type(points_file)==str
    point_or_index = 'OutputPoint'
    #get points
    with open(points_file, "r") as f:                
        lines=f.readlines()
        f.close()
    #####populate post-transformed array of contour centers
    sys.stdout.write("\n\n{} points detected\n\n".format(len(lines)))
    arr=np.empty((len(lines), 3))    
    for i in range(len(lines)):        
        arr[i,...]=lines[i].split()[lines[i].split().index(point_or_index)+3:lines[i].split().index(point_or_index)+6] #x,y,z            
    #optional save out of points
    dst_fl = os.path.join(dst, "posttransformed_zyx_voxels.npy")
    np.save(dst_fl, np.asarray([(z,y,x) for x,y,z in arr]))    
    #check to see if any points where found
    print("output array shape {}".format(arr.shape))
        
    return dst_fl

#%%
#path to points
pnts = arr
#transform
#make into transformix-friendly text file
transformed_dst = "/home/kepecs/Documents/cadaverine_slices/transform"
if not os.path.exists(transformed_dst): os.mkdir(transformed_dst)
pretransform_text_file = create_text_file_for_elastix(pnts, transformed_dst)
transformfiles = ["/home/kepecs/Documents/cadaverine_slices/elastix_inverse_transform/TransformParameters.0.txt",
                  "/home/kepecs/Documents/cadaverine_slices/elastix_inverse_transform/TransformParameters.1.txt"]
#copy over elastix files
transformfiles = modify_transform_files(transformfiles, transformed_dst) 
change_transform_parameter_initial_transform(transformfiles[0], 'NoInitialTransform')
#run transformix on points
points_file = point_transformix(pretransform_text_file, transformfiles[-1], transformed_dst)
#convert registered points into structure counts
converted_points = unpack_pnts(points_file, transformed_dst)

#%%
#map counts to structure
conv_pnts = np.load("/home/kepecs/Documents/cadaverine_slices/transform/posttransformed_zyx_voxels.npy")
from allensdk.api.queries.ontologies_api import OntologiesApi
oapi = OntologiesApi()
structure_graph = oapi.get_structures_with_sets([1])  # 1 is the id of the adult mouse structure graph
from allensdk.core.structure_tree import StructureTree
from allensdk.core.reference_space_cache import ReferenceSpaceCache

reference_space_key = 'annotation/ccf_2017'
resolution = 25
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
# This removes some unused fields returned by the query
structure_graph = StructureTree.clean_structures(structure_graph)  
tree = StructureTree(structure_graph)
name_map = tree.get_name_map()
# get the ids of all the structure sets in the tree
structure_set_ids = tree.get_structure_sets()
rsp = rspc.get_reference_space()

# query the API for information on those structure sets
allen_stuff = pd.DataFrame(oapi.get_structure_sets(structure_set_ids))
target_group = 167587189 #summary structures
summary_structures = tree.get_structures_by_set_id([target_group])
annotation, meta = rspc.get_annotation_volume()
annsag = np.fliplr(np.rot90(np.swapaxes(annotation, 0,2), axes=(1,2))) #reorient to sagittal
#
conv_pnts = conv_pnts.astype(int)
#get ids of + points of vasculature
def iid(annotation, pnt):
    z,y,x = pnt
    try:
        return annotation[z,y,x]
    except Exception as e:
        print(e)
    
iids = [iid(annsag,pnt) for pnt in conv_pnts]
#vasculature df
iidsdf = pd.DataFrame(iids).value_counts() #pandas series

ontology_file = "/home/kepecs/Documents/allen.json"
def get_progeny(dic,parent_structure,progeny_list):
    if "msg" in list(dic.keys()): dic = dic["msg"][0]
    name = dic.get("id")
    children = dic.get("children")
    if name == parent_structure:
        for child in children: # child is a dict
            child_name = child.get("id")
            progeny_list.append(child_name)
            get_progeny(child,parent_structure=child_name,progeny_list=progeny_list)
        return
    
    for child in children:
        child_name = child.get("id")
        get_progeny(child,parent_structure=parent_structure,progeny_list=progeny_list)
    return 

with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

iidf = pd.DataFrame()
iidf["id"] = pd.DataFrame(summary_structures)["id"]
iidf["name"] = pd.DataFrame(summary_structures)["name"]
iidf["acronym"] = pd.DataFrame(summary_structures)["acronym"]

for soi in iidf["id"] :
    progeny = []; get_progeny(ontology_dict, soi, progeny)
    #try except statements are to make sure all sois and children are included
    try:
        counts = [iidsdf[soi]]              #look at vasculature df
    except:
        counts = []
    for prog in progeny:
        try: 
            counts.append(iidsdf[prog])
        except:
            counts.append(0)
    try: #if an entry already exists
        iidf.loc[iidf["id" ]== soi, "vasculature_points"] = iidf.loc[iidf["id"] == soi, "vasculature_points"].sum() + np.sum(counts)
    except:
        iidf.loc[iidf.id == soi, "vasculature_points"] = np.sum(counts)

iidf.to_csv("/home/kepecs/Documents/cadaverine_slices/vasculature.csv", index=None)

#%%
#check
conv_pnts_ = conv_pnts[np.where((conv_pnts[:,0] > 100) & (conv_pnts[:,0] < 300))]
mapann = np.zeros_like(annsag)
for pnt_ in conv_pnts_:
    z,y,x = pnt_
    try:
        mapann[z,y,x] = 10000
    except Exception as e:
        print(e)
plt.imshow(np.max(mapann[100:300], axis = 0))    