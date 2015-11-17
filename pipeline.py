from __future__ import print_function
import os
import datetime
import ftplib
from scipy import *
import pylab as plt
import numpy as np
import healpy as hp
import readTEC

import pylab, os, sys
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

# Generate the necessary IONEX filename for a given day, and fetch it via ftp
def IONEX_file_needed(year, month, day):
	time_str = '{year} {month} {day}'.format(year=year, month=month, day=day)
	day_of_year = datetime.datetime.strptime(time_str, '%Y %m %d').timetuple().tm_yday

	if day_of_year < 10:
		day_of_year = '00{day_of_year}'.format(day_of_year=day_of_year)
	elif 10 <= day_of_year < 100:
		day_of_year = '0{day_of_year}'.format(day_of_year=day_of_year)

	# Outputing the name of the IONEX file you require
	IONEX_file = 'CODG{day_of_year}0.{year_end}I'.format(day_of_year=day_of_year, year_end=str(year)[2:4])

	return IONEX_file

def get_IONEX_file(year, month, day):
	server = 'ftp.unibe.ch'

	ftp_dir = os.path.join('aiub/CODE/', year)
	IONEX_file = IONEX_file_needed(year, month, day)
	IONEX_fileZ = ''.join((IONEXfile, '.Z'))

	getting_file_str = 'Retrieving {IONEX_fileZ} for {day} {month} {year}'.format(IONEX_fileZ=IONEX_fileZ, day=day, month=month, year=year)
	print(getting_file_str)

	ftp = ftplib.FTP(server, 'anonymous', 'jaguirre@sas.upenn.edu')
	ftp.cwd(ftp_dir)
	ftp.retrbinary(' '.join(('RETR', IONEX_fileZ)), open(IONEX_fileZ, 'wb').write)
	ftp.quit()

	return True

#------------------------------------------
# This script reads a IONEX file and retrieves
# the height of the Ionosphere
#------------------------------------------

def calc_ion_height(filename): 
	# opening and reading the IONEX file into memory
	with open(filename, 'r') as read_file:
		linestring = read_file.read()
		IONEX_list = linestring.split('\n')

	for file_data in IONEX_list:
		if file_data.split()[-1] == 'DHGT'
			ion_h = float(file_data.split()[0])

	return ion_h

# This function provides de geographic and topographic coordinates
# the Ionospheric piercing point (IPP)
def punc_ion_offset(lat_obs, az_sou, zen_sou, alt_ion):
	radius_earth = 6371000.0 # in meters

	# The 2-D sine rule gives the zenith angle at the
	# Ionospheric piercing point
	zen_punc = math.asin((radius_earth * math.sin(zen_sou)) / (radius_earth + alt_ion)) 

	# Use the sum of the internal angles of a triange to determine theta
	theta = zen_sou - zen_punc

	# The cosine rule for spherical triangles gives us the latitude
	# at the IPP
	lat_ion = math.asin(math.sin(lat_obs) * math.cos(theta) + math.cos(lat_obs) * math.sin(theta) * math.cos(az_sou)) 
	d_lat = lat_ion - lat_obs # latitude difference

	# Longitude difference using the 3-D sine rule (or for spherical triangles)
	d_lon = math.asin(math.sin(az_sou) * math.sin(theta) / math.cos(lat_ion))

	# Azimuth at the IPP using the 3-D sine rule
	s_az_ion = math.sin(az_sou) * math.cos(lat_obs) / math.cos(lat_ion)
	az_punc = math.asin(s_az_ion)

	return d_lat, d_lon, az_punc, zen_punc

def gen_hp_map(filename):
	filename = 'IONEX_Data/CODG3400.11I'

	test = readTEC.calc_TEC(-30.6988641207, 22.0338381509, filename)

	# In the current version, this gives just the 71 x 73 lat/lon array at
	# 13 values of UT.  calc_TEC performs an additional interpolation in
	# time to give values every hour.  There is also code to interpolates more finely spatially
	a = readTEC.read_IONEX_TEC(filename)

	# Now I have a standard problem solved before: grid a rectilinear
	# function of theta,phi onto a healpix map.  Ugh, why would IONEX pick
	# such an awful discretization?  It's oversampled at the poles, and sparse elsewhere
	# Index is lat, then lon

	# Original Healpix gridding
	nside = 32
	npix = hp.nside2npix(nside)

	grid_map = np.zeros(npix)
	hits = np.zeros(npix)

	# just straight-up HEALpix-ellize that bitch
	nside = 32
	npix = hp.nside2npix(nside)
	ipix = np.arange(npix)
	theta, phi = hp.pix2ang(nside, ipix)

	pix = hp.ang2pix(nside, theta, phi)
	#print 'Pix values'
	#print pix.min()
	#print pix.max()

	# Simplest gridding is
	#map[pix] = val
	# This tries to do some averaging
	for i, v in enumerate(val):
		grid_map[pix[i]] += v
		hits[pix[i]] += 1
	grid_map = grid_map / hits

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

def readmap():
	GENERATE = False
	TEST = False

	#2012-02-13T00:00:00 4h30m00s-30d43m17.5s CODG0440.12I.txt

	lat = -30.76528 * u.deg
	lon = 21.42831 * u.deg
	height = 1000 * u.m

	karoo = EarthLocation(lat=lat, lon=lon, height=height)
	date_str = '2012-02-13T22:00:00'
	time = Time('2012-02-13 22:00:00')
	PAPER_lat = '30d43m17.5ss'
	PAPER_lon = '21d25m41.9se'
	IONEX = 'CODG0440.12I'

	nside = 32
	npix = hp.nside2npix(nside)
	ipix = np.arange(npix)
	theta, phi = hp.pix2ang(nside, ipix)

	alt = (90. - np.degrees(np.array(theta))) * u.degree
	az = (np.degrees(np.array(phi))) * u.degree

	altaz = SkyCoord(alt=alt, az=az, obstime=time, frame='altaz', location=karoo)

	ra_h, ra_m, ra_s = altaz.icrs.ra.hms
	dec_h, dec_m, dec_s = altaz.icrs.dec.dms

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
		dec_str_base = '{dec_d}h{dec_m}m{dec_s}s'.format(dec_d=int(dec_d[p]), dec_m=int(abs(dec_m[p])), dec_s=abs(dec_s[p]))

		if dec_d[p] > 0:
			dec_str = ''.join(('+', dec_str_base))
		else:
			dec_str = dec_str_base

		if GENERATE:
			gen_str = 'IONFRM.py {ra_str}{dec_str} {PAPER_lat} {PAPER_lon} {date_str} {IONEX}'.format(ra_str=ra_str, dec_str=dec_str,
																									PAPER_lat=PAPER_lat, PAPER_lon=PAPER_lon,
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

if __name__ == '__main__':
	gen_hp_map(filename)
	readmap()
