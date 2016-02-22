from __future__ import print_function
import os
import sys
import numpy as np
import pylab as plt
import healpy as hp
import rad
from astropy import units as u
import time

def std_hour(UT):
    #print(int(UT))
    if UT < 10:
        hour = '0{hour}'.format(hour=int(UT))
    else:
        hour = '{hour}'.format(hour=int(UT))

    return hour

base_path = os.path.expanduser('~/radionopy')
#TEC, RMS, info = rad.read_IONEX_TEC('CODG1400.04I')

#toplot = []
for num in range(22, 23):
    #coord_file = os.path.join(base_path, 'RM_files/coords{num}.txt'.format(num=std_hour(num)))
    #coord_lat, coord_lon = np.loadtxt(coord_file, unpack=True)
    #nlat = len(coord_lat)
    #nlon = len(coord_lon)
    #lat_rad = np.outer(np.radians(90. - coord_lat), np.ones(nlon))
    #lon_rad = np.outer(np.ones(nlat), np.radians(coord_lon % 360))
    my_rad = os.path.join(base_path, 'RM_files/IonRM{num}.txt'.format(num=std_hour(num)))
    #my_rad = os.path.join(base_path, 'correct/IonRM{num}.txt'.format(num=std_hour(num)))
    UT, TEC, B, RM, dRM = np.loadtxt(my_rad, unpack=True)
    #checkmap = hp.get_interp_val(RM, lat_rad, lon_rad)
    print(RM.shape)
    #print(checkmap.shape)
    #toplot.append(checkmap)
    hp.orthview(RM, rot=[0, 90], min=0, max=2, title='SAST 00:00 2012-02-13', half_sky=True)
    #hp.mollview(RM, flip='geo')
plt.show()

#for i, t in enumerate(toplot):
#    plt.figure(i)
#    plt.clf()
#    plt.imshow(t)
#    plt.colorbar()

#plt.show()
