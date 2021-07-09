# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 11:18:50 2021

@author: zahhr
"""

import os, pandas as pd, matplotlib.pyplot as plt, seaborn as sns, numpy as np, tifffile as tif

src = r"C:\Users\zahhr\Box\kepecs_lab_summer2021\cachexia_blood_brain_barrier\images_staining_pilot\quantification"
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

#images to quantify
df = pd.DataFrame()
imgs = [os.path.join(src,xx) for xx in os.listdir(src) if "40x" in xx and "npy" not in xx and "threshold" not in xx]    
print(imgs)
for img in imgs:
    im = tif.imread(img)
    #0=dapi, 1=cldn5, 2=igg
    imgt = find_site(im[2], thresh=2, filter_kernel=(3,3)) #threshold by IgG (to find vasculature?) instead
    plt.imshow(imgt, "gist_yarg");plt.axis("off")
    plt.savefig(os.path.join(src, "threshold_"+os.path.basename(img)))
    # np.save(os.path.join(src, os.path.basename(img)+"_IgG_binary_threshold_IgG.npy"), imgt)
    #apply threshold to org image GREEN CHANNEL ONLY
    thres = im[2][imgt]
    #calculate stats
    df.loc[os.path.basename(img),"total thresholded pixels"] = imgt.sum() #to normalize
    df.loc[os.path.basename(img),"integrated intensity"] = thres.sum()
    df.loc[os.path.basename(img),"mean intensity"] = thres.mean()
    df.loc[os.path.basename(img),"std intensity"] = thres.std()
    df.loc[os.path.basename(img),"median intensity"] = np.median(thres)
#make group and structure labels
df["label"] = np.ones((len(df)))*np.nan
df.loc[df.index.str.contains("nh"), "label"] = "control"
df.loc[~df.index.str.contains("nh"), "label"] = "cachexia"
df["structure"] = np.ones((len(df)))*np.nan
df.loc[df.index.str.contains("vta"), "structure"] = "vta"
df.loc[df.index.str.contains("nac"), "structure"] = "nac"
df.loc[df.index.str.contains("hippocampus"), "structure"] = "hippocampus"
df.to_csv(os.path.join(src, "igg_igg_mask_python_threshold_analysis.csv"))
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
