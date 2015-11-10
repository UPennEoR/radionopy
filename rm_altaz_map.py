import numpy as np
import pylab as plt
import healpy as hp
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

karoo = EarthLocation(lat=-30.76528*u.deg,lon=21.42831*u.deg,height=1000*u.m)
time = Time('2011-09-15 04:44:42.481956')

nside=32
npix = hp.nside2npix(nside)
ipix = np.arange(npix)
theta,phi = hp.pix2ang(nside,ipix)

alt = (90.-np.degrees(np.array(theta)))*u.degree
az = (np.degrees(np.array(phi)))*u.degree

altaz = SkyCoord(alt=alt,az=az,obstime=time,frame='altaz',location=karoo)




