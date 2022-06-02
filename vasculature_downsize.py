#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 12:22:29 2022

@author: kepecs
"""

import os, numpy as np, tifffile as tif, SimpleITK as sitk, cv2, multiprocessing as mp, shutil, sys
from scipy.ndimage import zoom


def resize_helper(im, i, pth, dst, resizef):
    y,x = im.shape
    yr = int(y/resizef); xr = int(x/resizef)
    im = cv2.resize(im, (xr, yr), interpolation=cv2.INTER_LINEAR)
    tif.imsave(os.path.join(dst, os.path.basename(pth)[:-4]+str(i).zfill(4)+".tif"), 
                    im.astype("uint8"), compress=1)

if __name__ == "__main__":
    
    pth = "/home/kepecs/BALBc-no1_iso3um_stitched_segmentation.tif"
    print("\nPath to stitched images: %s\n" % pth)
    #path to store downsized images
    dst = "/home/kepecs/Documents/cadaverine_slices/vasculature_reg"
    print("\nPath to storage directory: %s\n\n" % dst)
    if not os.path.exists(dst): os.mkdir(dst)
    img = tif.memmap(pth)
    z = img.shape[0]
    resizef = 5 #factor to downsize imgs by
    iterlst = [(img[ii], ii, pth, dst, resizef) for ii in range(len(img))]
    p = mp.Pool(12)
    p.starmap(resize_helper, iterlst)
    p.terminate()
    
    #now downsample to 140% of allen atlas
    imgs = [os.path.join(dst, xx) for xx in os.listdir(dst) if "tif" in xx]; imgs.sort()
    z = len(imgs)
    y,x = sitk.GetArrayFromImage(sitk.ReadImage(imgs[0])).shape
    arr = np.zeros((z,y,x))
    #get allen atlas
    atlpth = "/home/kepecs/Documents/cadaverine_slices/average_template_25.tif"
    atl = sitk.GetArrayFromImage(sitk.ReadImage(atlpth))
    atlz,atly,atlx = atl.shape #get shape, sagittal
    #read all the downsized images
    for i,img in enumerate(imgs):
        if i%5000==0: print(i)
        arr[i,:,:] = sitk.GetArrayFromImage(sitk.ReadImage(img)) #horizontal
    #switch to sagittal
    arrsag = np.swapaxes(arr,2,0)
    z,y,x = arrsag.shape
    print((z,y,x))
    print("\n**********downsizing....heavy!**********\n")
    
    arrsagd = zoom(arrsag, ((atlz*1.4/z),(atly*1.4/y),(atlx*1.4/x)), order=1)
    tif.imsave(os.path.join(os.path.dirname(dst), "downsized_for_atlas.tif"), arrsagd.astype("uint16"))
    print("\ndeleting storage directory after making volume...\n %s" % dst)
    shutil.rmtree(dst)
