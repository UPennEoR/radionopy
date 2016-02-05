# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 09:08:04 2016

@author: jaguirre
"""

import numpy as np
import pylab as plt
import healpy as hp
import healpyTools as hpt
import polarized_beams as pb
reload(pb)

nside = 32
npix = hp.nside2npix(nside)
ipix = np.arange(npix)
t,p = hp.pix2ang(nside,ipix)

# Theta hat components in Cartesian basis
thx = np.cos(t)*np.cos(p)
thy = np.cos(t)*np.sin(p)
thz = -np.sin(t)
theta_hat_cart = np.stack((thx,thy,thz))

# Phi hat components in Cartesian basis
phx = -np.sin(p)
phy = np.cos(p)
phz = np.zeros_like(p)
phi_hat_cart = np.stack((phx,phy,phz))

# Define the rotation
rot = [0,90]
R = hp.Rotator(rot=rot)

# Perform the rotation
theta_hat_rot_cart = R(theta_hat_cart)
phi_hat_rot_cart = R(phi_hat_cart)

# Extract the components of the rotated vector in the original theta_hat,phi_hat
# basis
R_tt = np.einsum('a...,a...',theta_hat_rot_cart,theta_hat_cart)
R_tp = np.einsum('a...,a...',theta_hat_rot_cart,phi_hat_cart)
R_pt = np.einsum('a...,a...',phi_hat_rot_cart,theta_hat_cart)
R_pp = np.einsum('a...,a...',phi_hat_rot_cart,phi_hat_cart)

# OK, now ... might it do what I want?
zy = pb.zy_ideal_dipole(nside=32)
xy = pb.xy_ideal_dipole(nside=32)
pb.plot_efield_beams(zy['xt'],zy['xp'],zy['yt'],zy['yp'],figno=1)
pb.plot_efield_beams(xy['xt'],xy['xp'],xy['yt'],xy['yp'],figno=2)

zyr = {}
for k in zy.keys():
    zyr[k] = hpt.rotate_healpix_map(zy[k],[0,-90])

zy_rxt = zyr['xt'] * R_tt + zyr['xp'] * R_tp 
zy_rxp = zyr['xt'] * R_pt + zyr['xp'] * R_pp 
zy_ryt = zyr['yt'] * R_tt + zyr['yp'] * R_tp 
zy_ryp = zyr['yt'] * R_pt + zyr['yp'] * R_pp 

pb.plot_efield_beams(zy_rxt,zy_rxp,zy_ryt,zy_ryp,figno=3)



