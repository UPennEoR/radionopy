from __future__ import print_function
import os
import sys
import numpy as np
import healpy as hp
import pylab
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

def orth(base_UT22, start_UT22, unit, show=False):
    UT22_interp = np.zeros_like(base_UT22)
    c = 0
    for i, val in enumerate(start_UT22):
        if np.isnan(val):
            c += 1
            theta, phi = hp.pix2ang(nside, i)
            neybs = hp.get_neighbours(nside, theta, phi=phi) 
            UT22_interp[i] = np.nanmean(start_UT22[neybs[0]])

        else:
            UT22_interp[i] = val

    print(c, 'NaN vals')

    hp.orthview(UT22_interp, rot=[0, 90], min=0, max=2, unit=unit, title='SAST 00:00 2012-02-13', half_sky=True)
    if show:
        pylab.show()

    return UT22_interp

if __name__ == '__main__':
    GENERATE = False
    TEST = False

    #2012-02-13T00:00:00 4h30m00s-30d43m17.5s CODG0440.12I.txt

    lat = -30.76528 * u.deg
    lon = 21.42831 * u.deg
    height = 1000 * u.m

    karoo = EarthLocation(lat=lat, lon=lon, height=height)
    date_str = '2012-02-13T22:00:00'
    time = Time('2012-02-13 22:00:00')
    PAPERlat = '30d43m17.5ss'
    PAPERlon = '21d25m41.9se'
    IONEX = 'CODG0440.12I'

    nside = 32
    npix = hp.nside2npix(nside)
    ipix = np.arange(npix)
    theta, phi = hp.pix2ang(nside, ipix)

    alt = (90. - np.degrees(np.array(theta))) * u.degree
    az = (np.degrees(np.array(phi))) * u.degree

    altaz = SkyCoord(alt=alt, az=az, obstime=time, frame='altaz', location=karoo)

    ra_h, ra_m, ra_s = altaz.icrs.ra.hms
    dec_d, dec_m, dec_s = altaz.icrs.dec.dms

    storage = np.zeros((npix, 15, 2)) #number of pixels, 15 UTs each, RMs and eRMs for each UT

    UT0 = np.zeros((npix))
    UT4 = np.zeros((npix))
    UT6 = np.zeros((npix))
    UT22 = np.zeros((npix))
    eUT22 = np.zeros((npix))

    '''
    UT1_11pm
    UT2_12am
    UT3_1am
    UT4_2am
    '''

    c = 0
    for p in range(npix):
        ra_str = '{ra_h}h{ra_m}m{ra_s}s'.format(ra_h=int(ra_h[p]), ra_m=int(ra_m[p]), ra_s=int(ra_s[p]))
        dec_str_base = '{dec_d}d{dec_m}m{dec_s}s'.format(dec_d=int(dec_d[p]), dec_m=int(abs(dec_m[p])), dec_s=abs(dec_s[p]))

        if dec_d[p] > 0:
            dec_str = ''.join(('+', dec_str_base))
        else:
            dec_str = dec_str_base

        if GENERATE:
            gen_str = 'IONFRM.py {ra_str}{dec_str} {PAPERlat} {PAPERlon} {date_str} {IONEX}'.format(ra_str=ra_str, dec_str=dec_str,
                                                                                                    PAPERlat=PAPERlat, PAPERlon=PAPERlon,
                                                                                                    date_str=date_str, IONEX=IONEX)
            if TEST: 
                print(gen_str)
            else:
                os.system(gen_str)
        else:
            #READ
            filename = ''.join(('2012-02-13', ra_str, dec_str, 'IonRM.txt'))
            #print('Reading {filename}'.format(filename=filename))
            try: 
                UT, _, _, RM, eRM = np.loadtxt(filename, unpack=True)
                #print(p, UT.shape)
                try:
                    for i, ut in enumerate(UT):
                        if ut == 0.:
                            UT0[p] = RM[i]
                        elif ut == 4.:
                            UT4[p] = RM[i]
                        elif ut == 6.:
                            UT6[p] = RM[i]
                        elif ut == 22.: 
                            UT22[p] = RM[i]
                            eUT22[p] = eRM[i]
                except TypeError:
                    continue
            except IOError:
                print('issue with {filename}'.format(filename=filename))
                UT22[p] = np.nan
                c += 1
                continue

    print(c, 'IOErrors')
    hp.orthview(UT22, rot=[0, 90], max=2, unit=r'rad m$^{-1}$', title='SAST 00:00 2012-02-13', half_sky=True)
    #pylab.show()

    print(UT22)
    print('Interpolating NaNs')

    UT22_interp = orth(UT22, UT22, unit=r'rad m$^{-1}$', show=False)
    UT22_interp_2 = orth(UT22, UT22_interp, unit=r'rad m$^{-2}$', show=True)

    np.savez('mapdata.npz', map=UT22_interp_2)
    hp.write_map('mapdata.fits', UT22_interp_2)
