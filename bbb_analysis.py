# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 11:18:50 2021

@author: zahhr
"""

import os, pandas as pd, matplotlib.pyplot as plt, seaborn as sns, numpy as np, tifffile as tif
import cv2
from skimage.morphology import ball

src = "/home/kepecs/Documents/slices/tifs/cadaverine/"
dst = "/home/kepecs/Documents/slices/tifs/cadaverine"

def find_site(im, thresh=10, filter_kernel=(5,5,5)):
    """Find a connected area of high intensity, using a basic filter + threshold + connected components approach
    
    by: bdeverett
    Parameters
    ----------
    img : np.ndarray
        3D stack in which to find site (technically need not be 3D, so long as filter parameter is adjusted accordingly)
    thresh: float
        threshold for site-of-interest intensity, in number of standard deviations above the mean
    filter_kernel: tuple
        kernel for filtering of image before thresholding    
    Returns
    --------
    bool array of volume where coordinates where detected
    """
    from scipy.ndimage.filters import gaussian_filter as gfilt
    from scipy.ndimage import label
    if type(im) == str: im = tif.imread(im)

    filtered = gfilt(im, filter_kernel)
    thresholded = filtered > filtered.mean() + thresh*filtered.std() 
    labelled,nlab = label(thresholded)

    sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
    vals = [i+1 for i in np.argsort(sizes)[::-1]]
    return np.in1d(labelled, vals).reshape(labelled.shape)

def perimeter_sphericity_intensity(src, dims = 3):
    """
    src = 3d
    
    looks at it from two perspectives and then takes the total average
    dims=(2,3) number of dimensions to look at
   
    sometime two contours are found on a zplane after labels - in this case take min, but could take average?
    """
    #initialise
    sphericities = []; perimeters = []; ydim = []; xdim = []; intensity = [];
    if isinstance(src, str):
        im = tif.imread(os.path.join(os.path.dirname(src),os.path.basename(src)[10:]))
        src = tif.imread(src)        
    if 0 in src.shape: 
        return 0, 0 #projects against empty labels    
    else:
        contours = findContours(src.astype("uint8"))
        for contour in contours:
            x, y, _, _ = cv2.boundingRect(contour)
            circ = circularity(contour)
            perimeter = cv2.arcLength(contour, True)
            intensity_ = im[1,y,x]
            sphericities.append(circ); perimeters.append(perimeter); xdim.append(x); ydim.append(y); intensity.append(intensity_)
            
        #return values
        return perimeters, sphericities, xdim, ydim, intensity


def circularity(c):
    """
    A Hu moment invariant as a shape circularity measure, Zunic et al, 2010
    """
    circ = (4*np.pi*cv2.contourArea(c))/(cv2.arcLength(c,True)**2)

    return circ

def findContours(z):
    
    contours,hierarchy = cv2.findContours(z, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE) #if need more than two values to unpack error here upgrade to cv2 3+   
    contours = np.asarray([c.squeeze() for c in contours if cv2.contourArea(c)>0])
    
    return contours
#%%
#images to quantify
# df = pd.DataFrame()
imgs = ["/home/kepecs/Documents/slices/tifs/cadaverine/20um_cadaverine_40x_quantification_cadw2rh_20um_series3_slide1_bregma074_nac_40x.tif", 
        "/home/kepecs/Documents/slices/tifs/cadaverine/20um_cadaverine_40x_quantification_pr2w2lh_20um_series1_slide_bregma134-15_nac_40x.tif", 
        "/home/kepecs/Documents/slices/tifs/cadaverine/20um_cadaverine_40x_quantification_pr2w2lh_20um_series1_slide1_bregma118_nac_40x.tif", 
        "/home/kepecs/Documents/slices/tifs/cadaverine/20um_cadaverine_40x_quantification_pr2w2rh_20um_series6_slide1_bregma134-15_nac_40x.tif", 
        "/home/kepecs/Documents/slices/tifs/cadaverine/20um_cadaverine_40x_quantification_pr2w2lh_20um_series1_slide1_hippocampus_40x.tif", 
        "/home/kepecs/Documents/slices/tifs/cadaverine/20um_cadaverine_40x_quantification_pr2w2lh_20um_series1_slide1_thalamus_40x.tif", 
        "/home/kepecs/Documents/slices/tifs/cadaverine/20um_cadaverine_40x_quantification_pr2w2lh_20um_series1_slide1_bregma098_nac_40x.tif", 
        "/home/kepecs/Documents/slices/tifs/cadaverine/20um_cadaverine_40x_quantification_pr2w2rh_20um_series6_slide1_bregma118_nac_40x.tif"]    
print(imgs)
for img in imgs:
    im = tif.imread(os.path.join(src,img))
    #0=dapi, 1=cldn5, 2=igg
    imgt = find_site(im[1], thresh=2, filter_kernel=(5,5)) #threshold by IgG (to find vasculature?) instead
    tif.imsave(os.path.join(dst, "threshold_"+os.path.basename(img)),imgt.astype("uint16"))
    # np.save(os.path.join(src, os.path.basename(img)+"_IgG_binary_threshold_IgG.npy"), imgt)
    #apply threshold to org image GREEN CHANNEL ONLY
    # thres = im[2][imgt]
    #calculate stats
    # df.loc[os.path.basename(img),"total thresholded pixels"] = imgt.sum() #to normalize
    # df.loc[os.path.basename(img),"integrated intensity"] = thres.sum()
    # df.loc[os.path.basename(img),"mean intensity"] = thres.mean()
    # df.loc[os.path.basename(img),"std intensity"] = thres.std()
    # df.loc[os.path.basename(img),"median intensity"] = np.median(thres)
#make group and structure labels
# df["label"] = np.ones((len(df)))*np.nan
# df.loc[df.index.str.contains("nh"), "label"] = "control"
# df.loc[~df.index.str.contains("nh"), "label"] = "cachexia"
# df["structure"] = np.ones((len(df)))*np.nan
# df.loc[df.index.str.contains("vta"), "structure"] = "vta"
# df.loc[df.index.str.contains("nac"), "structure"] = "nac"
# df.loc[df.index.str.contains("hippocampus"), "structure"] = "hippocampus"
# df.to_csv(os.path.join(src, "igg_igg_mask_python_threshold_analysis.csv"))

#get perimter, sphericity, and xydim of each contour
thres=['/home/kepecs/Documents/slices/tifs/cadaverine/threshold_cadaverine647_20um_sections_2021823_imp_ctrl_perfused2021730_cd647_nac_tile.tif', 
       '/home/kepecs/Documents/slices/tifs/cadaverine/threshold_cadaverine647_20um_sections_2021823_pr2w2lhrh_cd647_nac_tile.tif', 
       '/home/kepecs/Documents/slices/tifs/cadaverine/threshold_cadaverine647_20um_sections_2021823_pr2w2lhrh_cd647_nac_tile_2.tif', 
       '/home/kepecs/Documents/slices/tifs/cadaverine/threshold_cadaverine647_20um_sections_2021823_pr2w2rh_ctrl_cd647_nac_tile.tif'
       ]
for im in thres:
    df=pd.DataFrame()
    perimeter, sphericity, xdim, ydim, intensity = perimeter_sphericity_intensity(im,dims=2)
    df["x"] = xdim
    df["y"] = ydim
    df["perimeter"] = perimeter
    df["sphericity"] = sphericity
    df["intensity"] = intensity
    df.to_csv(os.path.join(dst, os.path.basename(im)[:-4]+".csv"), index = None)

#%%
#analyze/threshold
dfs = ["/home/kepecs/Documents/slices/tifs/cadaverine/threshold_cadaverine647_20um_sections_2021823_imp_ctrl_perfused2021730_cd647_nac_tile.csv",
       "/home/kepecs/Documents/slices/tifs/cadaverine/threshold_cadaverine647_20um_sections_2021823_pr2w2lhrh_cd647_nac_tile.csv", 
       "/home/kepecs/Documents/slices/tifs/cadaverine/threshold_cadaverine647_20um_sections_2021823_pr2w2lhrh_cd647_nac_tile_2.csv", 
       "/home/kepecs/Documents/slices/tifs/cadaverine/threshold_cadaverine647_20um_sections_2021823_pr2w2rh_ctrl_cd647_nac_tile.csv"]
datadf = pd.DataFrame()
datadf["name"] = [os.path.basename(dfpth)[46:-4] for dfpth in dfs]
datadf["cells"] = [0]*len(dfs)
for dfpth in dfs:
    df = pd.read_csv(dfpth, index_col = None)
    print(dfpth)
    print(df.shape)
    plt.figure()
    plt.hist(df.sphericity.values)
    plt.ylabel("sphericity")
    plt.figure()
    plt.hist(df.perimeter.values)
    plt.ylabel("perimeter")
    plt.figure()
    plt.hist(df.intensity.values)
    plt.ylabel("intensity")
    #filter by sphericity - remove rounded objects
    df = df[df.sphericity < 0.8]
    df = df[(df.perimeter > 25) & (df.perimeter < 60)] #get non-specly objects
    df = df[(df.intensity > 400)] #filter by intensity
    df.to_csv(dfpth[:-4]+"_filtered_sphericity0d8_perimeter20_intensity400.csv", index = None)
    cellmap = np.zeros_like(tif.imread(dfpth[:-4]+".tif"))
    for x,y in np.array([df["x"].values,df["y"].values]).T:
        cellmap[y,x] = 255
    #apply x y dilation
    r = 6
    selem = ball(r)[int(r/2)]
    cellmap = cv2.dilate(cellmap, selem, iterations = 1)
    tif.imsave(dfpth[:-4]+"_filtered_sphericity0d8_perimeter20_intensity400.tif", cellmap.astype("uint8"))
    datadf.loc[datadf["name"] == os.path.basename(dfpth)[46:-4], "cells"] = len(df)
#%%
#normalize
df["normalized integrated intensity"] = df["integrated intensity"]/df["total thresholded pixels"]
#plot?
plt.figure()
sns.stripplot(data=df, x="label", y="integrated intensity", hue="structure",size=8)
plt.figure()
sns.stripplot(data=df, x="label", y="mean intensity", hue="structure",size=8)
plt.figure()
sns.stripplot(data=df, x="label", y="median intensity", hue="structure",size=8)
#ttest for nac
from scipy.stats import ttest_ind
tstat,pval = ttest_ind(df.loc[(df.label=="control") & (df.structure=="nac"), "normalized integrated intensity"].values,
                       df.loc[(df.label=="cachexia") & (df.structure=="nac"), "normalized integrated intensity"].values)
print("pvalue: %0.10f" % pval)
