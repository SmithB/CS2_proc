#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 16:21:12 2023

@author: ben
"""

import numpy as np
import pointCollection as pc
import pyproj
from utils import gaussian


class CS2_L1(object):
    def __init__(self, P=None, C=None, Ph=None, noise_DB=None, sigma_k=4, hemisphere=1):
        self.P=P
        self.C=C
        self.Ph=Ph
        self.N_samps=None
        self.P_noise=None
        if P is not None:
            self.N_samps=P.shape[1]
            self.N_bursts=P.shape[0]
            self.P_noise=10**(noise_DB/10)
        self.hemisphere=hemisphere
        self.sigma_k=sigma_k
        self.K_smooth = gaussian(np.arange(-4*sigma_k, 4.01*sigma_k), 0, sigma_k)
        self.K_smooth_phase = gaussian(np.arange(-16, 16.01), 0, 5)
        self.K_smooth_slope = np.convolve(np.array([1, 0, -1])/2, self.K_smooth, mode='full')
        self.six_sigma_pad = np.arange(-np.floor(6*sigma_k), np.floor(6*sigma_k), dtype=int)
        self.Ph_edge_correction = np.convolve(np.ones_like(self.six_sigma_pad, dtype=float), self.K_smooth, mode='same')
        self.data={}
        self.filename=None
        
        # constants:
        self.lambda_c = 0.022084  #(m)
        self.baseline = 1.1676    #(m)
        self.h_orbit = 732.e3  #(m, apehelion)
        # nominal conversion from phase to off-nadir distance
        self.dYdPhi = self.h_orbit * self.lambda_c / 2 / np.pi / self.baseline

        
        # build kernel for parabolic refinement of index location:
        self.parabolic_delta_ind=np.array([-1., 0., 1.])
        G=np.c_[np.ones_like(self.parabolic_delta_ind), self.parabolic_delta_ind, self.parabolic_delta_ind**2]
        self.parabolic_Ginv=np.linalg.inv(G.T @ G) @ G.T
        self.parabolic_delta_ind=self.parabolic_delta_ind.astype(int)
        self.llh2xyz = pyproj.Transformer.from_crs(4326, 4978)
        if self.hemisphere==1:
            self.EPSG=3413
        elif hemisphere==-1:
            self.EPSG=3031
        # spacecraft position and basis vectors
        self.xyz=None
        
        # EPSG 4326 is latitude, longitude
        xyz_EPSG = 4978   # WGS84 xyz
        self.xforms = dict(
            xyh2xyz = pyproj.Transformer.from_crs(self.EPSG, xyz_EPSG),
            xyz2xyh = pyproj.Transformer.from_crs(xyz_EPSG, self.EPSG),
            xy2ll = pyproj.Transformer.from_crs(self.EPSG, 4326),
            ll2xy = pyproj.Transformer.from_crs(4326, self.EPSG),
            lle2xyz = pyproj.Transformer.from_crs(4326, xyz_EPSG),
        )

    def xyh2xyz(self, xyh):
        return np.c_[self.xforms['xyh2xyz'].transform(*xyh.T)]

    def xyz2xyh(self, xyz):
        return np.c_[self.xforms['xyz2xyh'].transform(*xyz.T)]
    
    def lle2xyz(self, lle):
        return np.c_[self.xforms['lle2xyz'].transform(*lle.T)]
    
    def from_nc(self, filename, sigma_k=4):
        
        fields_1hz=['iono_cor_01',
         'load_tide_01',
         'mod_dry_tropo_cor_01',
         'mod_wet_tropo_cor_01',
         'ocean_tide_eq_01',
         'pole_tide_01',
         'solid_earth_tide_01']
                    
        fields_20hz = ['alt_20_ku',
         'ind_meas_1hz_20_ku',
         'lat_20_ku',
         'lon_20_ku',
         'flag_mcd_20_ku',
         'noise_power_20_ku',
         'off_nadir_roll_angle_str_20_ku',
         'sat_vel_vec_20_ku',
         'rec_count_20_ku',
         'time_20_ku',
         'transmit_pwr_20_ku',
         'uso_cor_20_ku',
         'window_del_20_ku']
        
        fields_WF = [ 'ph_diff_waveform_20_ku', 'coherence_waveform_20_ku', 'pwr_waveform_20_ku']
        fields_scale = ['echo_scale_factor_20_ku','echo_scale_pwr_20_ku']
        
        D_1hz = pc.data().from_h5(filename, fields=fields_1hz)
        D_20hz = pc.data().from_h5(filename, fields=fields_20hz)
        D_WF = pc.data().from_h5(filename, fields=fields_WF)
        D_scale = pc.data().from_h5(filename, fields=fields_scale)
        
        D_WF.coherence_waveform_20_ku[D_WF.coherence_waveform_20_ku>2]=np.NaN    
        D_pwr = pc.data().from_h5(filename, fields=['pwr_waveform_20_ku'])
        D_WF.assign(power = D_pwr.pwr_waveform_20_ku * D_scale.echo_scale_factor_20_ku[:, None] * (2.**D_scale.echo_scale_pwr_20_ku[:, None]))
        self=CS2_L1( P=D_WF.power, 
                 C=D_WF.coherence_waveform_20_ku,
                 Ph=D_WF.ph_diff_waveform_20_ku, 
                 noise_DB=D_20hz.noise_power_20_ku, 
                 sigma_k = sigma_k, 
                 hemisphere=self.hemisphere)
        self.data={'WF':D_WF,'1Hz':D_1hz,'20Hz':D_20hz}
        self.filename=filename
        return self

    def get_SC_coords(self):
        
        D_20hz = self.data['20Hz']
        
        lle=np.c_[D_20hz.lat_20_ku, D_20hz.lon_20_ku, D_20hz.alt_20_ku]
        x_SC = self.lle2xyz(lle)
        # a point 1 km above the SC:
        # nadir vector:
        bvX = x_SC-self.lle2xyz(np.c_[D_20hz.lat_20_ku, D_20hz.lon_20_ku, D_20hz.alt_20_ku+1000])
        bvX /= np.sqrt(np.sum(bvX**2, axis=1))[:, None]
    
        # satellite velocity vector in ITRF
        vel_hat = D_20hz.sat_vel_vec_20_ku
        vel_hat /= np.sqrt(np.sum(vel_hat**2, axis=1))[:, None]
    
        # Ellipsoidal velocity vector
        dotXV = np.sum( bvX * vel_hat, axis=1)
        bvY = vel_hat - (bvX * dotXV[:,None])
        bvY /= np.sqrt(np.sum(bvY**2, axis=1))[:, None]
    
        #cross product of Y vector and X vector to find Z vector
        bvZ = np.c_[bvY[:,1]*bvX[:,2] - bvY[:,2]*bvX[:,1], 
             bvY[:,2]*bvX[:,0] - bvY[:,0]*bvX[:,2],    
             bvY[:,0]*bvX[:,1] - bvY[:,1]*bvX[:,0]]
        # check that this is the same:
        bvZ = np.cross(bvY, bvX)
        
        xyh = self.xyz2xyh(x_SC)
        
        self.xyz=pc.data(columns=3).from_dict(dict(
            bvX=bvX,
            bvY=bvY,
            bvZ=bvZ,
            lle=lle,
            xyh_SC=xyh,
            vel_hat=vel_hat,
            x_SC=x_SC            
            ))
        
        return

