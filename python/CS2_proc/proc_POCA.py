#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 16:47:56 2023

@author: ben
"""

import numpy as np
from .utils import  find_maxima, conv_same
import pointCollection as pc

def pick_bursts(D_L1, bursts=None, C_min=0.75):

    if bursts is None:
        bursts=range(D_L1.N_bursts)

    D_sub=[]
    for burst in range(D_L1.N_bursts):
        D_sub += [ pick_burst(D_L1, burst, C_min=C_min) ]
    D_sub = [Di for Di in D_sub if Di is not None]

    return pc.data().from_list(D_sub)

def pick_burst(D_L1, burst, C_min=0.75, DEBUG=False):

    P=D_L1.P[burst,:]
    dPdt = np.convolve(P, D_L1.K_smooth_slope,  mode='same')
    P_smooth = np.convolve(P, D_L1.K_smooth,  mode='same')
    C_smooth = np.convolve(D_L1.C[burst,:], D_L1.K_smooth, mode='same')

    mask = (dPdt > 6*D_L1.P_noise[burst]/np.sqrt(D_L1.sigma_k)) & (D_L1.C[burst,:] > C_min)

    # eliminate small islands in the mask:
    mask=np.convolve(
        np.convolve(mask, np.ones(D_L1.sigma_k+1),mode='same')>D_L1.sigma_k/2,
            np.ones(D_L1.sigma_k), 'same')>0;
    Ph=D_L1.Ph[burst,:]

    maxima = find_maxima(dPdt, mask)
    if len(maxima)==0:
        return

    P_maxima = find_maxima(P_smooth)
    P_maxima = P_maxima[P_smooth[P_maxima] > 0.5*np.sqrt(np.mean(P_smooth**2))];

    maxima = maxima[dPdt[maxima] > 0.5*(P_smooth[maxima]/2/D_L1.sigma_k)]
    P_envelope = np.interp( maxima, np.concatenate([[0], P_maxima]), P_smooth[np.concatenate([[P_maxima[0]], P_maxima])])
    maxima = maxima[P_smooth[maxima] > 0.2*P_envelope]
    if len(maxima)==0:
        return

    # label is a variable that increments every time the mask goes 0->1, and is -1 where mask == 0
    label=np.concatenate([[0], np.cumsum(np.diff(mask)==1)])
    label[mask==0]=-1

    D_list=[]
    for i_max, this_max_ind in enumerate(maxima):
        this_D_out = pick_segment(D_L1, i_max, this_max_ind,
                                    P, Ph, dPdt, P_smooth, C_smooth, mask, label,
                                      DEBUG=DEBUG)
        if this_D_out is not None:
            this_D_out['burst']=burst+np.zeros_like(this_D_out['power'], dtype=int)
            D_list += [this_D_out]

    D_out=None
    if len(D_list):
        D_out={}
        fields=D_list[0].keys()
        for field in fields:
            D_out[field]=np.c_[[Di[field] for Di in D_list]]
        D_out=pc.data().from_dict(D_out)
    return D_out

def pick_segment(D_L1, i_max, max_ind, P, Ph, dPdt, P_smooth, C_smooth, mask, label, DEBUG=False):

    # find the valid region around the maximum
    this_label=label[max_ind]
    these_samps=np.flatnonzero(label==label[max_ind])
    i01 = [these_samps[0], these_samps[-1]]

    # set Smax_ind to the unrefined maximum
    Smax_ind = int(max_ind)

    # do a parabolic refinement if there are enough samples in this window
    if (Smax_ind > i01[0]) & (Smax_ind < i01[1]):
        m = D_L1.parabolic_Ginv @ dPdt[Smax_ind + D_L1.parabolic_delta_ind]
        Smax_ind += m[1]/2/m[2]

    # calculate the phase for the region around Smax_ind
    reg1 = int(np.floor(Smax_ind)) + D_L1.six_sigma_pad
    reg1 = reg1[(reg1>=0) & (reg1 < D_L1.N_samps)]
    Ph1 = np.unwrap( Ph[reg1] )

    # anchor the unwrapping on the phase of the first bin in the mask in this region:
    delta_phi_unwrap = Ph1 - Ph[reg1]
    first_good_samp = np.argmin(label[reg1]==this_label)
    Ph1 -= delta_phi_unwrap[first_good_samp]

    if len(Ph1) == len(D_L1.Ph_edge_correction):
        Ph1s = np.convolve(Ph1, D_L1.K_smooth, mode='same')/D_L1.Ph_edge_correction
    else:
        if len(Ph1) < 3*D_L1.sigma_k:
            return
        Ph1s = conv_same(Ph1, D_L1.K_smooth)/ \
            conv_same(np.ones_like(Ph1), D_L1.K_smooth)

    # note that we are not (yet) calculating the values N_fp_bins, P_int, or power_before_pick

    # interpolate the smoothed values at the sample value
    try:
        D_out=dict(
            phase = np.interp(Smax_ind, reg1, Ph1s),
            power = np.interp(Smax_ind, reg1, P_smooth[reg1]),
            dPdt_est = np.interp(Smax_ind, reg1,  dPdt[reg1]),
            coherence = np.interp(Smax_ind, reg1, C_smooth[reg1]),
            P_raw = np.interp(Smax_ind, reg1, P[reg1]),
            ret_count = i_max,
            samp = Smax_ind)
    except Exception:
        pass

    if DEBUG:
        D_out.update(dict(
            reg1=reg1,
            Ph1=Ph1,
            Ph1s=Ph1s,
            dPdt=dPdt,
            Ph=Ph,
            P=P,
            mask=mask,
            P_smooth=P_smooth,
            C_smooth=C_smooth,
            label=label,
            max_ind = max_ind,
            mask_samps=these_samps,
        ))

    return D_out
