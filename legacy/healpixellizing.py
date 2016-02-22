from __future__ import print_function
import healpy as hp
import os
import sys
import rad
import numpy as np
import pylab as plt

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
            neighbours,weights = hp.get_interp_weights(nside,theta[i],phi[i])
            # Add weighted values to map
            map[neighbours] += v*weights
            # Keep track of weights
            hits[neighbours] += weights
        map = map/hits
        wh_no_hits = np.where(hits == 0)
        print('pixels with no hits', wh_no_hits[0].shape)
        map[wh_no_hits[0]] = hp.UNSEEN
    else:    
        for i,v in enumerate(f):
            map[pix[i]] += v
            hits[pix[i]] +=1
        map = map/hits

    return map

TEC,RMS,info = rad.read_IONEX_TEC('CODG1400.04I')
nlat = len(TEC['lat'])
nlon = len(TEC['lon'])

lat_rad = np.outer(np.radians(90.-TEC['lat']),np.ones(nlon))
lon_rad = np.outer(np.ones(nlat),np.radians(TEC['lon']%360))
TECmap = TEC['TEC'][0,:,:]
RMSmap = RMS['TEC'][0,:,:]

nside = 16
map = healpixellize(TEC['TEC'][0,:,:],lat_rad,lon_rad,nside)
ipix = np.arange(hp.nside2npix(nside))
t,p = hp.pix2ang(nside,ipix)

wh = (np.where(map == np.nan))[0]
for i,w in enumerate(wh):
    #neighbors = hp.get_neighbours(nside,t[i],p[i])
    neighbors = hp.get_interp_weights(nside,t[i],p[i])
    map[w] = np.median(neighbors)

# We can write out the map
hp.write_map('TECmap.fits',map)

# We can also request a list of interpolated values anywhere, in particular at the original locations from the IONEX file
checkmap = hp.get_interp_val(map,lat_rad,lon_rad)

toplot = [TECmap,checkmap,TECmap-checkmap,np.log10(np.abs(100.*(TECmap-checkmap)/TECmap)),RMSmap]
toplot_names = ['TEC','interp(TEC)','TEC - interp(TEC)','log(|% diff.|)','RMS(TEC)']
for i,t in enumerate(toplot):
    plt.figure(i)
    plt.clf()
    plt.imshow(t)
    plt.suptitle(toplot_names[i])
    plt.colorbar()

plt.show()

#plt.imshow(TECmap)
#hp.mollview(map,flip='geo')
#plt.show()




