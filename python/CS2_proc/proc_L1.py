#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 13:00:43 2023

@author: ben
"""

import numpy as np
import pointCollection as pc
import pyproj
from utils import gaussian

from CS2_L1 import CS2_L1
from proc_POCA import pick_bursts
from proc_swath import proc_swath_bursts, reduce_swath
import sys
from locate_CS import locate_returns
from select_ambiguity import get_DEM, select_ambiguity
import re
import os

def __main__():
    
    from argparse import ArgumentParser
    
    parser=ArgumentParser('reduce CS2 L1 file to POCA and swath-based elevation estimates')
    
    parser.add_argument('--hemisphere','-H', type=int, default=-1)
    parser.add_argument('--DEM_file','-D', type=str, required=True)
    parser.add_argument('L1_file', type=str, required=True)
    parser.add_argument('--out_base','-b', type=str, default='.')
    args=parser.parse_args()
    
    #TBD: parse filename (e.g. CS_LTA__SIR_SIN_1B_20180710T101404_20180710T101529_D001.nc)
    
    # make output directories
    CS2_re=re.compile('SIN_1B_(\d{4})(\d\d)(\d\d)T')
    year=CS2_re.search(args.L1_file).group(1)
    
    
    dirs={name:os.join([args.out_base]+[subs]) for name, subs in zip(
        ['top','year','POCA','swath'], 
        [[],  [year], [year, 'POCA'], [year,'swath']])}           
    # make output directories
    for thedir in dirs.values():        
        if not os.path.isdir(thedir):
            os.mkdir(os.path.join([args.out_base] + thedir))
    
    out_files={'swath':os.path.join(dirs['swath'], os.path.basename(args.L1_file).replace('.nc','_swath.h5')),
                'POCA':os.path.join(dirs['POCA'],  os.path.basename(args.L1_file).replace('.nc','_POCA.h5'))}

    # write output as tiled H5s
    
    
    D_L1 = CS2_L1(hemisphere=args.hemisphere).from_nc(args.L1_file)
    # assign time (note that this is TAI, not UTC, different from UTC by +32 seconds)
    D_L1.data['20Hz'].assign(time=D_L1.data['20Hz'].time_20_ku/24/3600/365.25+2000)
    D_L1.get_SC_coords()
    
    D_ambig=[]
   
    D_sw = proc_swath_bursts(D_L1)
    D_ambig += [locate_returns(D_L1, D_sw)]
    
    D_POCA = pick_bursts(D_L1, C_min=0.75)
    D_ambig += [locate_returns(D_L1, D_POCA)]

    get_DEM(D_ambig, args.DEM_file)
    
    D_sw = select_ambiguity(D_sw, D_ambig[0])
    D_POCA = select_ambiguity(D_POCA, D_ambig[1])
    D_POCA.assign(time=D_L1.data['20Hz'][D_POCA.burst.astype(int)])
    pc.indexedH5(bin_W=[1.e5, 1.e5]).to_h5(D_POCA, out_files['POCA'])
    #D_L1.xyz.to_h5(out_file, group='xyz', replace=False)

    D_sw_r = reduce_swath(D_sw, D_L1.dYdPhi)
    D_sw_r.assign(time=D_L1.data['20Hz'][D_sw_r.burst.astype(int)])
    # write the output file
    pc.indexedH5(bin_W=[1.e5, 1.e5]).to_h5(D_sw_r, out_files['swath'])
    
if __name__=='__main__':
    __main__()