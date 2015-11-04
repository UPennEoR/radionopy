# ionFR is really badly written ... its nominal modules (like rdalaz)
# are all too wrapped up with the arguments passed to the script to
# really understand what it's looking for.

import pylab as plt

import numpy as np
import healpy as hp


# It looks like calcTEC is actually self-contained.  It does, however,
# implicitly assume you want each hour over a day, and just one
# coordinate point.  Moreover, every call reloads all the GIM maps for
# a day, so it's super inefficient.
import readTEC
reload(readTEC)

filename='IONEX_Data/CODG3400.11I'

test = readTEC.calcTEC(-30.6988641207, 22.0338381509, filename)

# In the current version, this gives just the 71 x 73 lat/lon array at
# 13 values of UT.  calcTEC performs an additional interpolation in
# time to give values every hour.  There is also code to interpolates more finely spatially
a = readTEC.readIonexTEC(filename)

# Now I have a standard problem solved before: grid a rectilinear
# function of theta,phi onto a healpix map.  Ugh, why would IONEX pick
# such an awful discretization?  It's oversampled at the poles, and sparse elsewhere
# Index is lat, then lon

# Original Healpix gridding
nside=32
npix = hp.nside2npix(nside)

map = np.zeros(npix)
hits = np.zeros(npix)

# just straight-up HEALpix-ellize that bitch
pix = hp.ang2pix(nside,theta,phi)
#print 'Pix values'
#print pix.min()
#print pix.max()

# Simplest gridding is
#map[pix] = val
# This tries to do some averaging
for i,v in enumerate(val):
    map[pix[i]] += v
    hits[pix[i]] +=1
map = map/hits
