import numpy as np
import healpy as hp
import os
import zm_tools as zt

import polarized_beams as pb
import healpyTools as hpt


def dbi2linear(dbi):
    linear = np.power(10., dbi/10.)
    return linear

def get_hera_jones():
    cstfileX = 'HERA_DISH_paper_feed_cyl36_150mhz_X.txt'
    cstfileY = 'HERA_DISH_paper_feed_cyl36_150mhz_Y.txt'
    
    dataX = np.loadtxt(cstfileX,skiprows=2)
    dataY = np.loadtxt(cstfileY,skiprows=2)
    
    theta = np.radians(dataX[:,0])
    phi = np.radians(dataY[:,1])
    
    x_dir = dbi2linear(dataX[:,2])
    x_dir_max = x_dir.max()
    x_dir /= x_dir_max
    
    y_dir = dbi2linear(dataY[:,2])
    y_dir_max = y_dir.max()
    y_dir /= y_dir_max
    
    x_theta_dir = dbi2linear(dataX[:,3])/x_dir_max
    x_phi_dir = dbi2linear(dataX[:,5])/x_dir_max
    x_theta_phase = np.radians(dataX[:,4])
    x_phi_phase = np.radians(dataX[:,6])
    
    y_theta_dir = dbi2linear(dataY[:,3])/y_dir_max
    y_phi_dir = dbi2linear(dataY[:,5])/y_dir_max
    y_theta_phase = np.radians(dataY[:,4])
    y_phi_phase = np.radians(dataY[:,6])
    
    j_xt = np.sqrt(x_theta_dir) * np.exp(1.j * x_theta_phase )
    j_xp = np.sqrt(x_phi_dir) * np.exp(1.j * x_phi_phase)
    j_yt = np.sqrt(y_theta_dir) * np.exp(1.j * y_theta_phase)
    j_yp = np.sqrt(y_phi_dir) * np.exp(1.j * y_phi_phase)
    
    jones = np.array([[j_xt,j_xp],[j_yt,j_yp]])
    
    return {'jones':jones,'theta':theta,'phi':phi}

def get_hera_jones_healpixellized(nside):
    cstfileX = 'HERA_DISH_paper_feed_cyl36_150mhz_X.txt'
    cstfileY = 'HERA_DISH_paper_feed_cyl36_150mhz_Y.txt'
    
    nside_interp = 32 # to be used for the initial lo-res healpix gridding.
    
    dataX = np.loadtxt(cstfileX,skiprows=2)
    dataY = np.loadtxt(cstfileY,skiprows=2)
    
    theta = np.radians(dataX[:,0])
    phi = np.radians(dataY[:,1])
    
    x_dir = dbi2linear(dataX[:,2])
    x_dir_max = x_dir.max()
    x_dir /= x_dir_max
    
    y_dir = dbi2linear(dataY[:,2])
    y_dir_max = y_dir.max()
    y_dir /= y_dir_max
    
    x_theta_dir = zt.hires_healpixellize(dbi2linear(dataX[:,3])/x_dir_max,theta,phi,nside_interp,nside)
    x_phi_dir = zt.hires_healpixellize(dbi2linear(dataX[:,5])/x_dir_max,theta,phi,nside_interp,nside)
    x_theta_phase = zt.hires_healpixellize(np.radians(dataX[:,4]),theta,phi,nside_interp,nside)
    x_phi_phase = zt.hires_healpixellize(np.radians(dataX[:,6]),theta,phi,nside_interp,nside)
        
    y_theta_dir = zt.hires_healpixellize(dbi2linear(dataY[:,3])/y_dir_max,theta,phi,nside_interp,nside)
    y_phi_dir = zt.hires_healpixellize(dbi2linear(dataY[:,5])/y_dir_max,theta,phi,nside_interp,nside)
    y_theta_phase = zt.hires_healpixellize(np.radians(dataY[:,4]),theta,phi,nside_interp,nside)
    y_phi_phase = zt.hires_healpixellize(np.radians(dataY[:,6]),theta,phi,nside_interp,nside)
    
    j_xt = np.sqrt(x_theta_dir) * np.exp(1.j * x_theta_phase )
    j_xp = np.sqrt(x_phi_dir) * np.exp(1.j * x_phi_phase)
    j_yt = np.sqrt(y_theta_dir) * np.exp(1.j * y_theta_phase)
    j_yp = np.sqrt(y_phi_dir) * np.exp(1.j * y_phi_phase)
    
    jones_data = np.array([[j_xt,j_xp],[j_yt,j_yp]])
    
    return {'jones':jones_data,'nside':nside}
    

def get_hera_jones_healpixellized_lores():
    cstfileX = 'HERA_DISH_paper_feed_cyl36_150mhz_X.txt'
    cstfileY = 'HERA_DISH_paper_feed_cyl36_150mhz_Y.txt'
    
    nside_interp = 64 # to be used for the initial lo-res healpix gridding.
    
    dataX = np.loadtxt(cstfileX,skiprows=2)
    dataY = np.loadtxt(cstfileY,skiprows=2)
    
    theta = np.radians(dataX[:,0])
    phi = np.radians(dataY[:,1])
    
    x_dir = dbi2linear(dataX[:,2])
    x_dir_max = x_dir.max()
    x_dir /= x_dir_max
    
    y_dir = dbi2linear(dataY[:,2])
    y_dir_max = y_dir.max()
    y_dir /= y_dir_max
    
    x_theta_dir = zt.healpixellize(dbi2linear(dataX[:,3])/x_dir_max,theta,phi,nside_interp)
    x_phi_dir = zt.healpixellize(dbi2linear(dataX[:,5])/x_dir_max,theta,phi,nside_interp)
    x_theta_phase = zt.healpixellize(np.radians(dataX[:,4]),theta,phi,nside_interp)
    x_phi_phase = zt.healpixellize(np.radians(dataX[:,6]),theta,phi,nside_interp)
        
    y_theta_dir = zt.healpixellize(dbi2linear(dataY[:,3])/y_dir_max,theta,phi,nside_interp)
    y_phi_dir = zt.healpixellize(dbi2linear(dataY[:,5])/y_dir_max,theta,phi,nside_interp)
    y_theta_phase = zt.healpixellize(np.radians(dataY[:,4]),theta,phi,nside_interp)
    y_phi_phase = zt.healpixellize(np.radians(dataY[:,6]),theta,phi,nside_interp)
    
    j_xt = np.sqrt(x_theta_dir) * np.exp(1.j * x_theta_phase )
    j_xp = np.sqrt(x_phi_dir) * np.exp(1.j * x_phi_phase)
    j_yt = np.sqrt(y_theta_dir) * np.exp(1.j * y_theta_phase)
    j_yp = np.sqrt(y_phi_dir) * np.exp(1.j * y_phi_phase)
    
    jones_data = np.array([[j_xt,j_xp],[j_yt,j_yp]])
    
    return {'jones':jones_data,'nside':nside_interp}

# This function produces the correct output which is a stack of Kronecker products
# from the stack of Jones matrices. But surely it can be implemented more efficiently
def loop_kron(jones):
    np.asarray(jones)
    length = jones.shape[-1]
    kjones = np.empty((4,4,length),dtype=np.complex128)
    for i in range(length):
        kjones[:,:,i] = np.kron(jones[:,:,i],np.conjugate(jones[:,:,i]))
        
    return kjones

def jones_to_mueller(jones):
    np.asarray(jones)
    kjones = loop_kron(jones)
    
    S = np.array([
            [1.,1.,0.,0.],
            [0.,0.,1.,1.j],
            [0.,0.,1.,-1.j],
            [1.,-1.,0.,0.]])/2.
    
    invS = np.array([
            [1.,0.,0.,1.],
            [1.,0.,0.,-1.],
            [0.,1.,1.,0.],
            [0.,-1.j,1.j,0.]])
    
    mueller = np.einsum('ab...,bc...,cd...->ad...',invS,kjones,S)
    return mueller

def jones_to_partial_mueller(jones):
    np.asarray(jones)
    kjones = loop_kron(jones)
    
    S = np.array([
            [1.,1.,0.,0.],
            [0.,0.,1.,1.j],
            [0.,0.,1.,-1.j],
            [1.,-1.,0.,0.]])/2.
    
    invS = np.array([
            [1.,0.,0.,1.],
            [1.,0.,0.,-1.],
            [0.,1.,1.,0.],
            [0.,-1.j,1.j,0.]])
    
    mueller = np.einsum('ab...,bc...->ac...',kjones,S)
    return mueller
