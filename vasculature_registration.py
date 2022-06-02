#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 12:59:43 2022

@author: kepecs
"""

import sys,os
sys.path.append("/home/kepecs/python/BrainPipe")
from tools.registration.register import elastix_command_line_call

src = "/home/kepecs/Documents/cadaverine_slices"
param_fld = "/home/kepecs/python/BrainPipe/parameterfolder" #change if using rat
atl = "/home/kepecs/Documents/cadaverine_slices/average_template_25.tif" #defaults to pma
    
#init registration
mv = os.path.join(src, "downsized_for_atlas.tif")
print("\nPath to downsized vol for registration to atlas: %s" % mv)
fx = atl
print("\nPath to atlas: %s" % fx)
out = os.path.join(os.path.dirname(src), "elastix")
if not os.path.exists(out): os.mkdir(out)

params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
#run
e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)


#atlas to registration vol
#inverse transform
fx = os.path.join(src, "downsized_for_atlas.tif")
mv = atl
print("\nPath to downsized vol for inverse registration to atlas: %s" % fx)
print("\nPath to atlas: %s" % mv)
out = os.path.join(src, "elastix_inverse_transform")
if not os.path.exists(out): os.mkdir(out)

params = [os.path.join(param_fld, xx) for xx in os.listdir(param_fld)]
#run
e_out, transformfiles = elastix_command_line_call(fx, mv, out, params)
