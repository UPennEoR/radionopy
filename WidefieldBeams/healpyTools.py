# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 14:40:50 2013

@author: jaguirre
"""

import healpy as hp
import numpy as np  

def rotate_healpix_map(map,rot):
    """ Will rotate the pixels of a map into (effectively) a new ordering representing a rotation of the function.  Not sure why this isn't implemented in healpy directly (maybe it is).  In order to map each pixel exactly to a new one, the transform is only accurate to the pixel size.  """
    
    npix = len(map)
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

def healpixellize(f_in,theta_in,phi_in,nside,fancy=True):
    """ A dumb method for converting data f sampled at points theta and phi (not on a healpix grid) into a healpix at resolution nside """

    # Input arrays are likely to be rectangular, but this is inconvenient
    f = f_in.flatten()
    theta = theta_in.flatten()
    phi = phi_in.flatten()
    
    pix = hp.ang2pix(nside,theta,phi)

    map = np.zeros(hp.nside2npix(nside))
    hits = np.zeros(hp.nside2npix(nside))
    
    # Simplest gridding is map[pix] = val. This tries to do some
    #averaging Better would be to do some weighting by distance from
    #pixel center or something ...
    if (fancy):
        for i,v in enumerate(f):
            # Find the nearest pixels to the pixel in question
            neighbours,weights = hp.get_neighbours(nside,theta[i],phi[i])
            # Add weighted values to map
            map[neighbours] += v*weights
            # Keep track of weights
            hits[neighbours] += weights
        map = map/hits
        wh_no_hits = np.where(hits == 0)
        print 'pixels with no hits',wh_no_hits[0].shape
        map[wh_no_hits[0]] = hp.UNSEEN
    else:    
        for i,v in enumerate(f):
            map[pix[i]] += v
            hits[pix[i]] +=1
        map = map/hits

    return map
