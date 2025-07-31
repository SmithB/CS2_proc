#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 13:00:43 2023

@author: ben
"""

import pointCollection as pc
import CS2_proc as CS2p
import sys
import re
import os
import glob

def make_filenames(args, L1_file):

    # Outputs will be stored by year
    CS2_re=re.compile('SIN_1B_(\d{4})(\d\d)(\d\d)T')
    year=CS2_re.search(L1_file).group(1)

    # make output directories
    dirs={}
    for name, subs in zip( ['top','year','POCA','swath'],
                          [[''],  [year], [year, 'POCA'], [year,'swath']]):
        temp=[args.retrack_base]+subs
        dirs[name]=os.path.join(*temp)

    for thedir in dirs.values():
        if not os.path.isdir(thedir):
            os.mkdir(thedir)

    # define output files
    out_files={'swath':os.path.join(dirs['swath'], os.path.basename(L1_file).replace('.nc','_swath.h5')),
                'POCA':os.path.join(dirs['POCA'],  os.path.basename(L1_file).replace('.nc','_POCA.h5'))}
    return out_files

def make_queue(args, argv):

    base_cmd = ' '.join([arg for arg in argv if arg not in ['--queue_glob']+['-q'] +args.queue_glob ])

    for queue_str in args.queue_glob:
        for L1_file in glob.glob(queue_str):
            if not args.Replace:
                # check if output files exist
                out_files = make_filenames(args, L1_file)
                if os.path.isfile(out_files['POCA']) and os.path.isfile(out_files['swath']):
                    continue
            print(f"{base_cmd} --L1_file {L1_file} ")

def __main__():

    from argparse import ArgumentParser

    parser=ArgumentParser('reduce CS2 L1 file to POCA and swath-based elevation estimates',
                          fromfile_prefix_chars="@")
    parser.add_argument('--hemisphere','-H', type=int, default=-1)
    parser.add_argument('--DEM_file','-D', type=str, required=True)
    parser.add_argument('--L1_file', type=str)
    parser.add_argument('--queue_glob','-q', type=str, nargs='+')
    parser.add_argument('--retrack_base','-b', type=str, default='.')
    parser.add_argument('--Replace', action='store_true')
    args, _=parser.parse_known_args()

    if args.queue_glob is not None:
        make_queue(args, sys.argv)
        return


    out_files = make_filenames(args, args.L1_file)
    for out_file in out_files.values():
        if os.path.isfile(out_file):
            os.remove(out_file)

    # read the L1 file
    D_L1 = CS2p.CS2_L1(hemisphere=args.hemisphere).from_nc(args.L1_file)
    # assign time (note that this is TAI, not UTC, different from UTC by +32 seconds)
    D_L1.data['20Hz'].assign(time=D_L1.data['20Hz'].time_20_ku/24/3600/365.25+2000)
    D_L1.get_SC_coords()

    # process the swath data, store in D_ambig
    D_ambig=[]
    D_sw = CS2p.proc_swath_bursts(D_L1)
    D_ambig += [CS2p.locate_returns(D_L1, D_sw)]

    # process the POCA data, store in D_ambig
    D_POCA = CS2p.pick_bursts(D_L1, C_min=0.75)
    D_ambig += [CS2p.locate_returns(D_L1, D_POCA)]

    # read the DEM for the swath and POCA data
    CS2p.get_DEM(D_ambig, args.DEM_file)

    # select ambiguities based on the DEM
    D_sw = CS2p.select_ambiguity(D_sw, D_ambig[0])
    D_POCA = CS2p.select_ambiguity(D_POCA, D_ambig[1])

    # write the POCA data to an indexed h5 (100-km bins)
    if D_POCA is None:
        print(f"No POCA points found for {args.L1_file}")
    else:
        D_POCA.assign(time=D_L1.data['20Hz'].time[D_POCA.burst.astype(int)])
        pc.indexedH5.data(bin_W=[1.e5, 1.e5]).to_file(D_POCA, out_files['POCA'])

    if D_sw is None:
        print(f"No swath points found for {args.L1_file}")
    else:
        # reduce the swath data to 1 point every ~200 m
        D_sw_r = CS2p.reduce_swath(D_sw, D_L1.dYdPhi)
        D_sw_r.assign(time=D_L1.data['20Hz'].time[D_sw_r.burst.astype(int)])
        # write the swath data to an indexed h5 (100-km bins)
        pc.indexedH5.data(bin_W=[1.e5, 1.e5]).to_file(D_sw_r, out_files['swath'])

if __name__=='__main__':
    __main__()
