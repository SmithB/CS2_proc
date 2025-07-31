#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 11:33:39 2024

@author: ben
"""
import numpy as np
import pointCollection as pc

def get_DEM(DAs, DEM_file):

    res=1.e4
    pad=5.e2
    ixys = []
    all_xy=set()
    for DA in DAs:
        _, ixy = pc.unique_by_rows( np.round(np.c_[DA.x[:,1], DA.y[:,1]]/res)*res, return_dict=True)
        ixys += [ixy]
        DA.assign(DEM=np.zeros_like(DA.x)+np.NaN)
        all_xy.update(ixy)

    for xy0 in all_xy:
        xy_temps = []
        to_search=[]
        for count, ixy in enumerate(ixys):
            if xy0 in ixy:
                xy_temps += [DAs[count].copy_subset(ixy[xy0], by_row=True, datasets=['x','y'])]
                to_search += [count]
        xy_temp = pc.data().from_list(xy_temps)
        DEM=pc.grid.data().from_geotif(DEM_file, bounds=xy_temp.bounds(pad=pad))

        for ind, xy_temp in zip(to_search, xy_temps):
            ii=ixys[ind][xy0]
            if hasattr(DEM,'z'):
                DAs[ind].DEM[ii,:]=DEM.interp(xy_temp.x, xy_temp.y)
            else:
                DAs[ind].DEM[ii,:]=np.NaN

def select_ambiguity(D0, DA):
    _, iBR = pc.unique_by_rows(np.c_[D0.burst, D0.ret_count], return_dict=True)
    r=DA.h-DA.DEM
    rows=[]
    cols=[]
    for BR, ind in iBR.items():
        ri = r[ind,:]
        keep = np.all(np.isfinite(ri), axis=1)
        if not np.all(keep):
            ri=ri[keep,:]
            ind=ind[keep]
        if len(ind)==0:
            continue

        R = np.mean(np.abs(ri), axis=0)
        if len(ind)==1:
            best = np.argmin(R)
        else:
            med_r = np.median(ri, axis=0)
            spread_r = np.mean(np.abs(ri-med_r), axis=0)

            best_R=np.argmin(R)
            best_spread=np.argmin(spread_r)
            if best_R==best_spread:
                best = best_R
            else:
                frac_R = (R-R.min())/np.maximum(2, R.min())
                frac_spread = (spread_r-spread_r.min())/np.maximum(1, spread_r.min())

                if frac_spread[best_R] < 0.5:
                    # We know that the best spread is worse than the spread for best R.  If the
                    # increment by which the spread for the best R is less than 50%, use the best_R
                    best = best_R
                elif frac_R[best_spread] < 0.1 :
                    # If the fractional disadvantage in R for the best spread value is less than 10%
                    # select the value associated with the best R
                    best = best_spread
                else:
                    # Both the spread and the R are significantly worse when the other determines the
                    # optimal ambiguity.  Don't use this segment.
                    continue
        rows += [ind]
        cols += [np.zeros_like(ind)+best]

    if len(rows)==0:
        return

    rows=np.concatenate(rows)
    cols=np.concatenate(cols)

    D_out=D0[rows].assign({field:getattr(DA, field)[rows, cols] for field in DA.fields})
    # add the phase shift associated with the best-fitting ambiguity to the phase
    delta_ph_ambig = 2*np.pi*np.array([-1, 0, 1])
    D_out.phase += delta_ph_ambig[cols]

    return D_out
