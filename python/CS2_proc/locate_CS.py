#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 13:16:10 2023

@author: ben
"""
import numpy as np
import pyproj
import pointCollection as pc

def locate_ambiguity(D_L1,  burst, range_surf, phase, roll):


    alpha = D_L1.lambda_c * phase / (2*np.pi*D_L1.baseline); #alpha = inferred angle (rad) 
 
    inf_angle = (alpha+roll)[:,None] #roll corrected inferred angle (rad)   !!! changed to + here!!!

    #calculate range vector
    Rhat = D_L1.xyz.bvX[burst,:] * np.cos(inf_angle) + D_L1.xyz.bvZ[burst,:] * np.sin(inf_angle)
       
    #new points in cartesian (m)
    x_new = D_L1.xyz.x_SC[burst,:] + Rhat*range_surf[:,None]
    
    # convert new points to lat/lon/h (rad,rad,m)
    xyh = D_L1.xyz2xyh(x_new)
  
    return xyh[:,0], xyh[:,1], xyh[:,2]



def locate_returns(D_L1b, D_ret, **params):

    D_1Hz=D_L1b.data['1Hz']
    D_20Hz=D_L1b.data['20Hz']
    D_WF=D_L1b.data['WF']

    if 'geod' not in params:
        params['geod']=pyproj.Geod(ellps='WGS84')

    D_20Hz.lat_20_ku[D_20Hz.lat_20_ku == 0] = np.NaN;
    D_20Hz.lon_20_ku[~np.isfinite(D_20Hz.lat_20_ku)] = np.NaN;
    
    #[CP.xps_sc,CP.yps_sc]=fwd_projection(D.lat_20_ku(:),D.lon_20_ku(:)); #polar stereographic
    #CP.h_sc=reshape(D.alt_20_ku(:), j, k);
    #calculate range
    c = 299792458; # speed of light (m/s)
    B = 3.2e8;     # measured chirp bandwidth (Hz)

    # NOTE: The USO correction is included in the window_del_20 field (for baselines D and E)
    # This was not apparent in the field description for L1D
    win_delay = D_20Hz.window_del_20_ku[D_ret.burst] #* (D_20Hz.uso_cor_20_ku[D_ret.burst]+1); # USO_Corr_factor, field 2
    N_samps=D_WF.power.shape[1]
    range_surf = (win_delay*c)/2 - (N_samps*c)/(8*B) + (D_ret.samp*c)/(4*B) #range(m)
    
    sec_num=D_20Hz.ind_meas_1hz_20_ku[D_ret.burst];
    
    for field in ['mod_dry_tropo_cor_01','mod_wet_tropo_cor_01','iono_cor_01','load_tide_01','solid_earth_tide_01','pole_tide_01']:
        range_surf += getattr(D_1Hz, field)[sec_num]

    D_ret.assign(range_surf = range_surf)
    
    # NOTE: baseline-D roll angle is in degrees
    roll=D_20Hz.off_nadir_roll_angle_str_20_ku[D_ret.burst]*np.pi/180
    delta_ph_ambig=2*np.pi*np.array([-1, 0, 1]);

    D_ambig = pc.data().from_dict({field:np.zeros((D_ret.samp.size, 3)) for field in ['x','y','h']})
    
    for col, dPhi in enumerate(delta_ph_ambig):
        D_ambig.x[:,col], D_ambig.y[:,col], D_ambig.h[:,col] = locate_ambiguity(D_L1b, D_ret.burst, D_ret.range_surf, D_ret.phase + dPhi, roll)

    
    return D_ambig

