import pylab, random, os, sys
import healpy as hp
import numpy as np
from astropy import units as u
from astropy import constants as c
from scipy import integrate
from bm_prms import prms

def rotate_hmap(map,rot):
	npix = map.shape[0]
	nside = hp.npix2nside(npix)

	rotmap = np.zeros(npix)
	ipix = np.arange(npix)
	t,p = hp.pix2ang(nside,ipix)

	r = hp.Rotator(rot=rot)

	# For each pixel in the new map, find where it would have come 
	# from in the old    
	trot,prot = r(t,p)
	ipix_rot = hp.ang2pix(nside,trot,prot)

	rotmap = map[ipix_rot]

	return rotmap

beam = np.load('XX_beam_maps.npz')['maps']
freqs = np.linspace(0.117,0.182,num=131) #aipy likes GHz units. avoiding band edges
nside = 128
npix = hp.nside2npix(nside)

beam_167 = beam[100,:]
#beam_167 += 1e-9


def brp(beam,rotarr,log=False):
	a1 = np.sqrt(beam)
	a2 = np.sqrt(rotate_hmap(beam_167,rotarr))
	if not log: return a1*a2
	else: return np.log10(a1*a2)

L = True
lb = '\n'
#f,axarr = pylab.subplots(4,1)

for i,rot in enumerate([0,10,20,30]):
	R = hp.orthview(brp(beam_167,[rot,0],log=L), rot=[180,0], half_sky=True, max=0, min=-3, title='Rotation=%i'%rot, unit=lb+r'$\log_{10}({\rm Normalized Response})$')#,hold=True)


pylab.show()