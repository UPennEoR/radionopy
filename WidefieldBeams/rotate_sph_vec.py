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

def vrot(vector,theta,phi,rot):
    """ 
        Given a 3-D unit vector specifying the vector field at (theta,phi),
        rotate this field by rot and return the new vector at the new point 
    """
    R = hp.Rotator(rot=rot)
    rvector = R(vector)
    rtheta,rphi = R(theta,phi)
    return (rvector,rtheta,rphi)

def t_hat_cart(t,p):
    """ Calculate the theta_hat vector at a given point (t,p) """
    thx = np.cos(t)*np.cos(p)
    thy = np.cos(t)*np.sin(p)
    thz = -np.sin(t)
    return np.stack((thx,thy,thz))

def p_hat_cart(t,p):
    """ Calculate the phi_hat vector at a given point (t,p) """
    phx = -np.sin(p)
    phy = np.cos(p)
    phz = np.zeros_like(p)
    return np.stack((phx,phy,phz))

def r_hat_cart(t,p):
    """ Calculate the r_hat vector at a given point (t,p) """
    rx = np.sin(t)*np.cos(p)
    ry = np.sin(t)*np.sin(p)
    rz = np.cos(t)
    return np.stack((rx,ry,rz))

def adot(x,y):
    """ The dot product term by term of a list of vectors """
    return np.einsum('a...,a...',x,y)

rot = [0,-90]

# Check a few things
t1 = 0. #np.pi/2.
p1 = 0. #np.pi/2.
th1 = t_hat_cart(t1,p1)
ph1 = p_hat_cart(t1,p1)
th1r,t1r,p1r = vrot(th1,t1,p1,rot)
ph1r,t1r,p1r = vrot(ph1,t1,p1,rot)
print np.round(adot(th1r,t_hat_cart(t1r,p1r)),decimals=2), np.round(adot(th1r,p_hat_cart(t1r,p1r)),decimals=2),np.degrees(t1r),np.degrees(p1r)
print np.round(adot(ph1r,t_hat_cart(t1r,p1r)),decimals=2), np.round(adot(ph1r,p_hat_cart(t1r,p1r)),decimals=2),np.degrees(t1r),np.degrees(p1r)

# Define Healpix coordinates
nside = 32
npix = hp.nside2npix(nside)
ipix = np.arange(npix)
t,p = hp.pix2ang(nside,ipix)

# Generate theta_hat and phi_hat in Cartesian basis
theta_hat_cart = t_hat_cart(t,p)
phi_hat_cart = p_hat_cart(t,p)

# Perform the rotation of the basis vectors
theta_hat_rot_cart,t_rot,p_rot = vrot(theta_hat_cart,t,p,rot)
phi_hat_rot_cart,t_rot,p_rot = vrot(phi_hat_cart,t,p,rot)
#t_rot,p_rot = R(t,p)

R_tt = adot(theta_hat_rot_cart,t_hat_cart(t_rot,p_rot))
R_tp = adot(theta_hat_rot_cart,p_hat_cart(t_rot,p_rot))
R_pt = adot(phi_hat_rot_cart,t_hat_cart(t_rot,p_rot))
R_pp = adot(phi_hat_rot_cart,p_hat_cart(t_rot,p_rot))

maps = [R_tt,R_tp,R_pt,R_pp]

if False:
    for map in maps:
        hp.orthview(map,rot=[0,90],half_sky=True)
    plt.show()

#
## OK, now ... might it do what I want?
zy = pb.zy_ideal_dipole(nside=32)
xy = pb.xy_ideal_dipole(nside=32)

pb.plot_efield_beams(zy['xt'],zy['xp'],zy['yt'],zy['yp'],figno=1)
pb.plot_efield_beams(xy['xt'],xy['xp'],xy['yt'],xy['yp'],figno=2,rot=[0,90])

zyr = {}
for k in zy.keys():
    zyr[k] = hpt.rotate_healpix_map(zy[k],rot)

zy_rxt = zyr['xt'] * R_tt + zyr['xp'] * R_tp 
zy_rxp = zyr['xt'] * R_pt + zyr['xp'] * R_pp 
zy_ryt = zyr['yt'] * R_tt + zyr['yp'] * R_tp 
zy_ryp = zyr['yt'] * R_pt + zyr['yp'] * R_pp 

pb.plot_efield_beams(zy_rxt,zy_rxp,zy_ryt,zy_ryp,figno=3,rot=[0,90])

plt.show()

