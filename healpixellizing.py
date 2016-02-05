import healpy as hp
#import os
#import sys
path='/Users/jaguirre/PyModules/radionopy/'
sys.path.append(""+str(path)+"WidefieldBeams")
import healpyTools as hpt
reload(hpt)
import rad
import numpy as np

TEC,RMS,info = rad.read_IONEX_TEC('CODG1400.04I')
nlat = len(TEC['lat'])
nlon = len(TEC['lon'])

lat_rad = np.outer(np.radians(90.-TEC['lat']),np.ones(nlon))
lon_rad = np.outer(np.ones(nlat),np.radians(TEC['lon']%360))
TECmap = TEC['TEC'][0,:,:]
RMSmap = RMS['TEC'][0,:,:]

nside = 16
map = hpt.healpixellize(TEC['TEC'][0,:,:],lat_rad,lon_rad,nside)
ipix = np.arange(hp.nside2npix(nside))
t,p = hp.pix2ang(nside,ipix)

wh = (np.where(map == np.nan))[0]
for i,w in enumerate(wh):
    neighbors = hp.get_neighbours(nside,t[i],p[i])
    map[w] = np.median(neighbors)

# We can write out the map
hp.write_map('TECmap.fits',map)

# We can also request a list of interpolated values anywhere, in particular at the original locations from the IONEX file
checkmap = hp.get_interp_val(map,lat_rad,lon_rad)

toplot = [TECmap,checkmap,TECmap-checkmap,(TECmap-checkmap)/TECmap,RMSmap]
for i,t in enumerate(toplot):
    plt.figure(i)
    plt.clf()
    plt.imshow(t)
    plt.colorbar()

plt.show()

#plt.imshow(TECmap)
#hp.mollview(map,flip='geo')
#plt.show()




