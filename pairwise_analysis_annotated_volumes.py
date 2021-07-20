

import os, numpy as np, pandas as pd, matplotlib.pyplot as plt, tifffile as tif
from tools.conv_net.utils.io import pairwise_distance_metrics, read_roi_zip

if __name__ == "__main__":

    #paths and imports
    src = "/home/kepecs/Documents/cfos_annotations/volumes"
    dst = "/home/kepecs/Documents/cfos_annotations/results"
    brains = [os.path.join(src, xx) for xx in os.listdir(src) if not "." in xx]; brains.sort()
    
    imgsrc = "/home/kepecs/Documents/cfos_annotations/volumes"
    vols = [os.path.join(imgsrc, xx) for xx in os.listdir(imgsrc) if "tif" in xx]; vols.sort()
    roipths = [os.path.join(imgsrc, xx) for xx in os.listdir(imgsrc) if "RoiSet.zip" in xx]; roipths.sort()
    #these will be zyx 
    #note that importing it this way, the z dimension does not start from 0
    annotated_cells = np.array([np.array([[int(yy) for yy in xx.replace(".roi", "").split("-")] for xx in
                                read_roi_zip(roipth, points=True)]) for roipth in roipths], dtype = object)
    
    #voxel pair cutoff
    cutoffs = [5] #if a cell is more than this distance away from another cell, it is a separate cell
    
    #make dest
    if not os.path.exists(dst): os.mkdir(dst)
    
    for cutoff in cutoffs:
        for br, brain in enumerate(brains):
            fls = [os.path.join(brain, xx) for xx in os.listdir(brain)]; fls.sort()
            #init dataframe
            print("brain: %s\ncutoff: %s \n\n" % (os.path.basename(brain), cutoff))
            df = pd.DataFrame()
            df["parameters"] = [os.path.basename(xx) for xx in fls]
            #init
            df["tp"] = np.zeros(len(df))
            df["fp"] = np.zeros(len(df))
            df["fn"] = np.zeros(len(df))
            df["f1"] = np.zeros(len(df))
            
            for n, fl in enumerate(fls):
                if n%50 == 0: print(n)
                detected_cells = np.load(fl, allow_pickle=True)
                #change cells to zyx
                detected_cells =  np.array([detected_cells[c] for c in 'zyx']).T
                #generate overlays to inspect
                imbr = tif.imread(vols[br])
                cellmap = np.zeros_like(imbr)
                for pnt in detected_cells:
                    cellmap[pnt[0],pnt[1],pnt[2]] = 1
                merged = np.stack([imbr, cellmap, np.zeros_like(cellmap)], -1)
                tif.imsave(fl[:-4]+".tif", merged.astype("uint16"))
                detected_cells[:,0] = detected_cells[:,0]+1 #add 1 for roi comparison to accurate
                paired, tp, fp, fn = pairwise_distance_metrics(annotated_cells[br], 
                            detected_cells, cutoff = cutoff, verbose = False) 
                
                try: #calculating precision and recall
                    precision = tp/(tp+fp)
                    recall = tp/(tp+fn) #same as true positve rate
                    f1 = 2*( (precision*recall)/(precision+recall) ) #calculating f1 score
                except Exception as e:
                    print(e)
                    f1 = np.nan #if tp, fn, etc. are 0                
                    
                df.loc[df.parameters == os.path.basename(fl), "f1"] = f1
                df.loc[df.parameters == os.path.basename(fl), "tp"] = tp
                df.loc[df.parameters == os.path.basename(fl), "fp"] = fp
                df.loc[df.parameters == os.path.basename(fl), "fn"] = fn
                df.loc[df.parameters == os.path.basename(fl), "precision"] = precision
                df.loc[df.parameters == os.path.basename(fl), "recall_tpr"] = recall 
        
            #export csv per brain/volume                
            df.to_csv(os.path.join(dst, "{0}_{1}.csv".format(os.path.basename(brain), cutoff)))