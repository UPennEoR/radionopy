import numpy as np
import aipy as a
import healpy as hp
from bm_prms import prms
from healpyTools import rotate_healpix_map as rotmap

#freqs = [0.150] #np.arange(0.115,0.190,0.005)#0.150
#nside = 256
#npix = hp.nside2npix(nside)

def aipy_hp_beam(freq,nside,y=False,efield=False):
    # currently only works with scalar frequencies
    #for freq in freqs:
    #bm = prms['beam'](np.array([freq]),nside=nside,lmax=20,mmax=20,deg=7)
    #bm.set_params(prms['bm_prms'])
    px = np.arange(hp.nside2npix(nside))
    theta,phi = hp.pix2ang(nside,px)
    B = aipy_beam(freq,theta,phi,efield=efield)
    if y:
        B = rotmap(B,[90,0])
    #xyz = hp.pix2vec(nside,px)
    #poly = np.array([h.map[px] for h in bm.hmap])
    #Axx = np.polyval(poly,freq)
    #Axx = np.where(xyz[-1] >= 0, Axx, 0)
    #Axx /= Axx.max()
    #Axx = Axx*Axx
    return B

def aipy_beam(freq,theta,phi,efield=False):
    # really only defined between 120 and 180 MHz
    # theta, phi in radians
    nside = 256
    bm = prms['beam'](np.array([freq]),nside=nside,lmax=20,mmax=20,deg=7)
    bm.set_params(prms['bm_prms'])
    px = hp.ang2pix(nside,theta,phi)
    xyz = hp.pix2vec(nside,px)
    poly = np.array([h.map[px] for h in bm.hmap])
    Axx = np.polyval(poly,freq)
    Axx = np.where(xyz[-1] >= 0, Axx, 0)
    Axx /= Axx.max()
    if (not efield):
        Axx = Axx*Axx
    return Axx
    


