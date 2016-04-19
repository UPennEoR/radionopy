'''
radionopy.orth

authors | James Aguirre, Immanuel Washington, Saul Kohn

purpose | script to plot RM data into healpix orthview
'''
from __future__ import print_function
import os
import sys
import time
import numpy as np
import pylab as plt
import healpy as hp
from astropy import units as u
import radiono as rad

if __name__ == '__main__':
    for num in range(22, 23):
        my_rad = os.path.join(rad.base_path, 'RM_files/IonRM{num}.txt'.format(num=rad.std_hour(num)))
        UT, TEC, B, RM, dRM = np.loadtxt(my_rad, unpack=True)
        print(RM.shape)
        hp.orthview(RM, rot=[0, 90], min=0, max=2, title='SAST 00:00 2012-02-13', half_sky=True)
    plt.show()
