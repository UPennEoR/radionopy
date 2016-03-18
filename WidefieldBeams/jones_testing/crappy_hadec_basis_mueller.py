import numpy as np
import healpy as hp
import zm_tools as zt
import instrumental_response_tensors as irt

jones_data = irt.get_hera_jones_healpixellized_lores()
jones = jones_data['jones']

nside = jones_data['nside'] # nside=64, no other option right now.
npix = hp.nside2npix(nside)
pixindex = np.arange(npix)
za,ha = hp.pix2ang(nside,pixindex) # radians
dec = np.deg2rad(90. - np.rad2deg(za)) # radians

lat0 = -26. #degrees
codec0 = 90. - lat0
R = hp.Rotator(rot=[0.,codec0,0.]) 
Rinv = R.get_inverse()
rotmat = np.asarray(R.mat)

theta,phi = Rinv(za,ha)
localidx = hp.ang2pix(nside,theta,phi)

basis_rot = zt.get_sphr_rotation_matrix(za,ha,rotmat)

window_idx = np.where(za < np.pi /2.) # treating za as the theta coordinate
beam_window = np.zeros(npix)
beam_window[window_idx] = 1.

jones_hd_ae = jones[:,:,localidx] 
beam_window0 = beam_window[localidx]

beam_window0_arr = np.repeat(beam_window0,4).reshape(npix,2,2).transpose((1,2,0))

jones_masked_hd_ae = jones_hd_ae * beam_window0_arr # response is zeroed outside the local horizon

jones_hd_hd = np.einsum('ab...,bc...->ac...',jones_masked_hd_ae,basis_rot)

##$$$$$$$$$$$$$  JUST LOOK DOWN HERE.
mueller_hd_hd = irt.jones_to_mueller(jones_hd_hd) # mueller_hd_hd.shape = (4,4,npix), nside=64

# - first dimension is row, second dimension is column, third is sky location
    #Ex. mueller_hd_hd[0,2,i] is the Q->I component at the i'th pixel.
# - mueller_hd_hd[:,:,i] is the Mueller matrix that operates on a vector (I,Q,U,V) defined in the Ha/Dec basis when the array is at a -26 degrees Latitude on the Earth.
# - These components can be rotated around the z-axis, no more basis changes required (preeeeeeety sure at this point)

np.savez('socalled_leakage_beams',mueller_hd_hd)
