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

def rotate_vector_field(f_theta,f_phi,theta,phi,rot):
    """ Given a vector field with theta-hat component f_theta and phi-hat 
    component f_phi, where each component has a value at each coordinate 
    (theta,phi), rotate this field according to the healpy rot specification
    and return the rotated field in terms of its components theta-hat, phi-hat
    at the same positions specified in (theta,phi) """
    # Generate theta_hat and phi_hat in Cartesian basis
    theta_hat_cart = t_hat_cart(theta,phi)
    phi_hat_cart = p_hat_cart(theta,phi)
    # Perform the rotation of the basis vectors.  This is slightly inefficient,
    # since the coordinates are rotated twice.
    theta_hat_rot_cart,t_rot,p_rot = vrot(theta_hat_cart,theta,phi,rot)
    phi_hat_rot_cart,t_rot,p_rot = vrot(phi_hat_cart,theta,phi,rot)
    # Dot the new and old systems
    R_tt = adot(theta_hat_rot_cart,t_hat_cart(t_rot,p_rot))
    R_tp = adot(theta_hat_rot_cart,p_hat_cart(t_rot,p_rot))
    R_pt = adot(phi_hat_rot_cart,t_hat_cart(t_rot,p_rot))
    R_pp = adot(phi_hat_rot_cart,p_hat_cart(t_rot,p_rot))
    # Rotate the scalar functions
    f_theta_r = hpt.rotate_healpix_map(f_theta,rot)
    f_phi_r = hpt.rotate_healpix_map(f_phi,rot)
    f_theta_new = R_tt * f_theta_r + R_tp * f_phi_r
    f_phi_new = R_pt * f_theta_r + R_pp * f_phi_r
    return (f_theta_new,f_phi_new)
    

