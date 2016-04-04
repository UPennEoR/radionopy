import numpy as np
import healpy as hp
import zm_tools as zt
import instrumental_response_tensors as irt

nside = 64
npix = hp.nside2npix(nside)
rdindx = np.arange(npix) # healpix pixel index for Ra/Dec coords
codec,ra = hp.pix2ang(nside,rdindx)
# codec = np.pi - dec
dec = np.pi - codec

# A constant sky in the za,ha basis and (ha,dec) coordinates

## A single point source at zenith:
if True:
    cele_pole_idx = [0,1,2,3]
    cp_theta,cp_phi = hp.pix2ang(nside,cele_pole_idx)
    z_theta,z_phi = R(cp_theta,cp_phi)
    z_idx = hp.ang2pix(nside,z_theta,z_phi)
    I = np.zeros(npix)
    I[z_idx] = 10.
    Q = 0. * np.ones(npix)
    U = 0. * np.ones(npix)
    V = 0. * np.ones(npix)
    
## A flat sky    
if False:
    I = 10. * np.ones(npix)
    Q = 0. * np.ones(npix)
    U = 0. * np.ones(npix)
    V = 0. * np.ones(npix)



# get an active rotation matrix to the z'-axis through (dec,ha)=(-26,0)
dec0 = -26. # degrees
ha0 = 0. # degrees
codec0 = 90. - dec0 # degrees

# fuuuuuuuuuuuuuuuck this
Rinv = hp.Rotator(rot=[0,-codec0,0])
R = R.get_inverse()
rotation_matrix = np.asarray(Rinv.mat)
"""
This gives the matrix-valued function that transforms a vector field in the Ha/Dec
basis to the Az/El basis at the location (dec0,ha0). It turns out this function is 
symmetric under rotations about the celestial Z-axis. This is due to the fact that 
the spherical unit vector fields have a cylindrical symmetry. Hence, the rotation
from the Ha/Dec basis to the Az/El basis is constant as a function of time.

The resulting arrays are indexed in Ra/Dec coordinates.
"""

sphr_basis_rotation = zt.get_sphr_rotation_matrix(codec,ra,rotation_matrix)

# get the HERA Jones matrix
jones_data = irt.get_hera_jones_healpixellized_lores() # healpix maps

jones_ae_ae = jones_data['jones'] # jones.shape = (2,2,npix). 

# We define these maps to be in local Az/El coordinates, so lets change that:

window_idx = np.where(codec < np.pi /2.) # treating codec as the theta coordinate
beam_window = np.zeros(npix)
beam_window[window_idx] = 1.

theta0,phi0 = R(codec,ra)
ae0index = hp.ang2pix(nside,theta0,phi0)

jones_hd_ae = jones_ae_ae[:,:,ae0index] 
beam_window0 = beam_window[ae0index]

beam_window0_arr = np.repeat(beam_window0,4).reshape(npix,2,2).transpose((1,2,0))

jones_masked_hd_ae = jones_hd_ae * beam_window0_arr


# This is the Jones matrix that maps from the sky in the Az/El basis to the beam-weighted sky in the instrumental basis at time t = ha+ra = 0. 
# ((This is probably not quite right. The point is that we're going to rotate the scalar components to simulate the rotation of the instrument under the sky.))
# However, the scalar components are to be considered functions of Ha/Dec.

# Maybe this is a class? Then it can have attributes like "Coordinate system" and "Vector basis" and methods to change between them, 'n stuff...

jones_hd_hd = np.einsum('ab...,bc...->ac...',jones_masked_hd_ae,sphr_basis_rotation)

sky_hd_hd = np.array([[I+Q,U -1.j *V],[U+1.j*V,I-Q]])

## Lets see about some propagators for one baseline

u,v = 50.,50.

s = np.array(hp.pix2vec(nside,np.arange(npix)))
b = np.outer(np.array([u,v]),np.ones(npix))

b_dot_s = np.einsum('a...,a...',s[:2,:],b)
nu = 150. * 10.**6. # Hz
c = 299792458. # m/s

propagator_ae = np.exp(1.j * nu * b_dot_s / c)
propagator_hd = propagator_ae[ae0index]
propagator_hd_arr = np.repeat(propagator_hd,4).reshape(npix,2,2).transpose((1,2,0))

holographic_sky_hd_instr = np.einsum('ab...,bc...,dc...->ad...',jones_hd_hd,sky_hd_hd,np.conjugate(jones_hd_hd))

vis_150_0 = np.sum(holographic_sky_hd_instr * propagator_hd_arr,axis=2)
