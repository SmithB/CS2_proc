#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 31 12:29:32 2023

@author: ben
"""

import numpy as np
import pointCollection as pc
from utils import med_ind_and_spread_for_segment
import scipy.sparse as sps
from hermite_poly_fit import hermite_design_matrix


def proc_swath_bursts(D_L1, bursts=None, Cmin=0.75):
    
    if bursts is None:
        bursts=np.arange(D_L1.data['WF'].shape[0], dtype=int)
    
    dYdPhi=D_L1.h_orbit*D_L1.lambda_c/2/np.pi/D_L1.baseline
    phase_gap_size=200/dYdPhi
    
    D_sw=[]    
    for burst in bursts:
        P = D_L1.P[burst,:]
        Ph = D_L1.Ph[burst,:]
        C = D_L1.C[burst,:].copy()
        temp=proc_swath_return(Ph, C, P, burst, D_L1.K_smooth_phase, dYdPhi, phase_gap_size)
        if temp is not None:
            D_sw += temp
    D_sw = pc.data().from_list(D_sw)
    return D_sw

def calc_smooth_phase(Ph, C, P, K):
    
    
    
    #eliminate low-signal high-coherence sections (where we suspect that part 
    # of the burst has no signal from one RX chain);
    Cbad=np.convolve( (C>0.995) | ((C>0.900) & (P<2e-17)), np.ones(7), mode='same')>0;
    C[Cbad]=0;
 
    Phc=np.exp(1j*Ph)
    Cw=np.convolve(C, K, mode='same')
    Phs=np.angle(np.convolve(Phc*C, K,mode='same')/Cw);
    Cs=np.convolve(C*C, K,mode='same')/Cw
    
    return Phs, Cs
        

def proc_swath_return(Ph0, C0, P, burst,  K,  dYdphi,  phase_gap_size, C_threshold=0.75, DEBUG=False):
    
    # make the smoothed phase for this return
    Phs, Cs = calc_smooth_phase(Ph0, C0, P, K)
    
    # find regions where C is greater than the threshold
    mask = np.convolve(Cs > C_threshold, np.ones(5), mode='same')>2;
    
    # unwrap the phase to look for discontinuities 
    temp = np.unwrap(Phs)
    temp[mask==0] = np.NaN
    
    # find segments with continuous phase
    dPhi = np.diff(temp)
    continuous_phase = np.zeros_like(Ph0, dtype=bool)
    continuous_phase[:-1] = np.isfinite(dPhi) & (np.abs(dPhi) < phase_gap_size)
    
    mask[continuous_phase==0]=0
    
    # label is a variable that increments every time the mask goes 0->1, and is -1 where mask == 0
    label=np.concatenate([[0], np.cumsum(np.diff(mask.astype(float))==1)])
    label[mask==0]=-1

    _, D_label = pc.unique_by_rows(label, return_dict=True)
    D_label.pop((-1,))
    seg_count=-1   
    
    max_ind = np.argmax(P)
    
    D=[];
    for val, ind in D_label.items():
        # skip small regions
        if len(ind) < 20 or val==-1:
            continue
        ind=np.sort(ind)

        # skip segments that contain the waveform maximum [NOT IN MATLAB]
        if max_ind >= ind[0] and max_ind <= ind[-1]:
            continue
        seg_count += 1
        # unwrap, starting at the 10th sample (avoids edge problems)
        Ph_seg = np.NaN+np.zeros_like(ind)
        Ph_seg[10:] = np.unwrap(Phs[ind[10:]]);
        Ph_seg[10::-1] = np.unwrap(Phs[ind[10::-1]]);
        delta_amb=Ph_seg-Phs[ind];
        D += [pc.data().from_dict(dict(
            burst = np.zeros_like(Ph_seg, dtype=int)+burst,
            samp = ind,
            phase = Ph_seg,
            phase_raw = Ph0[ind]+delta_amb,
            ret_count = seg_count + np.zeros_like(Ph_seg, dtype=int),
            coherence = Cs[ind],   
            coherence_raw = C0[ind],
            power = P[ind]
        ))]
    return D
    
def reduce_swath(D_sw, dYdPhi):
    
    _, i_BR = pc.unique_by_rows(np.c_[D_sw.burst, D_sw.ret_count], return_dict=True)
    
    iMed, spread, block_N = [], [], []
    
    for count, (BR, i_BR) in enumerate(i_BR.items()):
        if len(i_BR) < 2:
            continue
        i_BR = np.sort(i_BR)
        y_ph = dYdPhi*D_sw.phase[i_BR]
        
        YR=(np.min(y_ph)-50, np.max(y_ph))
        YW=np.ceil((YR[1]-YR[0])/400)*400+50
        y0= np.arange(YR[0], YR[0]+YW, 400)
    
        if len(y0)==2:
            G_h=np.c_[y0[1]-y_ph, y_ph-y0[0]]/(y0[1]-y0[0])
        else:
            G_h = hermite_design_matrix(y_ph, y0)
        m = sps.linalg.spsolve(G_h.T @ G_h, G_h.T @ D_sw.h[i_BR])
        h_est = G_h @ m
        r = D_sw.h[i_BR]-h_est
    
        for delta_x in [0, 200]:
            iMed0, spread0, N0 = med_ind_and_spread_for_segment(y_ph, r, 400, x0=delta_x)
            keep = N0>0
            if np.any(keep):
                iMed += [i_BR[iMed0[keep,:]]]
                spread += [spread0[keep]]
                block_N += [N0[keep]]
    
    iMed = np.concatenate(iMed, axis=0)
    spread = np.concatenate(spread, axis=0)
    block_N = np.concatenate(block_N, axis=0)
    D_sw_r = D_sw[iMed]
    for field in D_sw_r.fields:
        setattr(D_sw_r, field, np.mean(getattr(D_sw_r, field), axis=1))
    D_sw_r.__update_size_and_shape__()
    
    return D_sw_r

