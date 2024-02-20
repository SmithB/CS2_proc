#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 12:40:48 2023

@author: ben
"""


import numpy as np

def gaussian(x, c=0, sigma=1):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-1/2*((x-c)/sigma)**2)

def find_maxima(x, mask=None):

    if mask is not None:
        return 1 + np.flatnonzero((x[1:-2] > x[0:-3]) & (x[1:-2] > x[2:-1]) & mask[1:-2])
    else:
        return 1 + np.flatnonzero((x[1:-2] > x[0:-3]) & (x[1:-2] > x[2:-1]))

def bin_spread(z, ind, sorted=False):
    z0 = z[ind]
    if sorted:
        zs=z0
    else:
        ii=np.argsort(z0)
        inverse_ind=np.zeros_like(ii)
        inverse_ind[ii]=np.arange(len(ii), dtype=int)
        zs=z0[ii]
    N=len(z0)
    C=np.arange(0.5, N)/N
    pp=np.interp(np.array([0.14, 0.5, 0.86]), C, zs)
    # need to find the indexes bounding the median
    if np.mod(N,2)==0:
        iMed = np.floor(N/2).astype(int)+np.array([-1, 0])
    else:
        iMed = np.floor(N/2).astype(int)+np.array([0, 0])
    if sorted:
        return ind[iMed], (pp[2]-pp[0])
    else:
        return ind[ii[iMed]], (pp[2]-pp[0])

def med_ind_and_spread_for_segment(x, z, delta, x0=0):
    ii=np.argsort(x)
    inverse_ind=np.zeros_like(ii)
    inverse_ind[ii]=np.arange(len(x), dtype=int)
    xs=x[ii]
    zs=z[ii]
    x_bin = np.round((xs-x0)/delta)*delta+x0
    #print(f'x_bin={str(x_bin)}')
    ux_bin = np.unique(x_bin)
    # first and last indices in each bin
    ii_first = np.searchsorted(x_bin, ux_bin, side='left')
    ii_last  = np.searchsorted(x_bin, ux_bin, side='right')-1
    N=np.zeros(len(ux_bin), dtype=int)
    iMed = np.zeros((len(ux_bin), 2), dtype=int)
    spread = np.zeros(len(ux_bin), dtype=float)
    for count, (xi, iif, iil), in enumerate(zip(ux_bin, ii_first, ii_last)):
        these=np.arange(iif, iil+1, dtype=int)
        this_x = xs[these]
        #make sure that there are selected samples on both sides of the bin center:
        # Otherwise we will process bins overlapping only the start and end of the segment
        if not (np.any(this_x < xi) and np.any(this_x > xi)):
            continue
        N[count] = iil-iif+1
        iMed[count,:], spread[count] = bin_spread(zs, these, sorted=False)
        #print(f'\t\t x:{xs[np.arange(iif, iil+1, dtype=int)]}, \n\t\t z:{zs[np.arange(iif, iil+1, dtype=int)]}')
        #print([xi, iif, iil, xs[iMed[count,:]], zs[iMed[count,:]], iMed[count,:]])

    # convert the index on the sorted x_phase to the unsorted index
    iMed=ii[iMed]
    return iMed, spread, N

def conv_same(aa, bb):
    if len(bb) <= len(aa):
        return np.convolve(aa, bb, mode='same')
    cf=np.convolve(aa, bb, mode='full')
    Nbh=int(len(bb)/2)
    if np.mod(len(bb),2)==1:
        return cf[Nbh:-Nbh]
    else:
        return cf[Nbh-1:Nbh]
