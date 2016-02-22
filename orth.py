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
    if UT < 10:
        hour = '0{hour}'.format(hour=int(UT))
    else:
        hour = '{hour}'.format(hour=int(UT))

    return hour

base_path = os.path.expanduser('~/radionopy')

for num in range(22, 23):
    my_rad = os.path.join(base_path, 'RM_files/IonRM{num}.txt'.format(num=std_hour(num)))
    UT, TEC, B, RM, dRM = np.loadtxt(my_rad, unpack=True)
    print(RM.shape)
    hp.orthview(RM, rot=[0, 90], min=0, max=2, title='SAST 00:00 2012-02-13', half_sky=True)
plt.show()
