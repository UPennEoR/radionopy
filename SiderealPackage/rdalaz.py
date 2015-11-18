#!/usr/bin/env python

# Slight modified version of rdaa.py
#================================================================
# rdaa: Convert right ascension/declination to azimuth/altitude
# For documentation, see:
#	 http://www.nmt.edu/tcc/help/lang/python/examples/sidereal/ims/
#----------------------------------------------------------------
#================================================================
# Imports
#----------------------------------------------------------------
from __future__ import print_function
import sys
import re
import sidereal
from math import *
#================================================================
# Manifest consants
#----------------------------------------------------------------

SIGN_PAT = re.compile(r'[\-+]')
#----- main

def alaz(tim):
	'''
	Main program for rdaa.
	'''
	print('WHAT THE HOLY FUCK?')
	print(tim)
	#-- 1 --
	# [ if sys.argv contains a valid set of command line
	# arguments ->
	#	 ra_dec := the right ascension and declination as
	#				a sidereal.ra_dec instance
	#	 lat_lon := the observer's location as a
	#				 sidereal.lat_lon instance
	#	 dt := the observer's date and time as a
	#			 datetime.datetime instance
	# else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	print('How does this work?!')
	ra_dec, lat_lon, dt = check_args(tim)
	print('ra_dec lat_lon dt ', ra_dec, lat_lon, dt)
	#-- 2 --
	# [ if dt has no time zone information ->
	#	 utc := dt
	# else ->
	#	 utc := the UTC equivalent to dt ]
	if (dt.tzinfo is None) or (dt.utcoffset() is None):
		utc = dt
	else:
		utc = dt - dt.utcoffset()
	#-- 3 --
	# [ sys.stdout +:= local sidereal time for dt and lat_lon ]
	gst = sidereal.SiderealTime.from_datetime(utc)
	lst = gst.lst(lat_lon.lon)
	#############print 'Equatorial coordinates:', ra_dec
	#############print 'Observer's location:', lat_lon
	#############print 'Observer's time:', dt
	#############print 'Local sidereal time is', lst
	#-- 4 --
	# [ h := hour angle for ra_dec at time (utc) and longitude
	#		 (lat_lon.lon) ]
	h = ra_dec.hour_angle(utc, lat_lon.lon)

	#############print 'Hour Angle:', h*180.0/pi,'d'

	#-- 5 --
	# [ aa := horizon coordinates of ra_dec at hour angle h
	#		 as a sidereal.AltAz instance ]
	aa = ra_dec.alt_az(h, lat_lon.lat)

	#-- 6 --
	#############print 'Horizon coordinates:', aa
	# - - - c h e c k A r g s
	print('ra_dec', ra_dec)
	
	# changing latitude of the observer in degrees to radians
	lat_deg = lat_lon.__str__().split()[0].split('d')[0].split('[')[1]
	lat_min = lat_lon.__str__().split()[1].split(''')[0]
	lat_sec = lat_lon.__str__().split()[2].split(''')[0]
	lat_degs = float(lat_deg) + float(lat_min) / 60.0 + float(lat_sec) / 3600.0
	lat_rads = lat_degs * pi / 180.0
	# changing longitude of the observer in degrees to radians
	lon_deg = lat_lon.__str__().split()[5].split('d')[0]
	lon_min = lat_lon.__str__().split()[6].split(''')[0]
	lon_sec = lat_lon.__str__().split()[7].split(''')[0]
	lon_degs = float(lon_deg) + float(lon_min) / 60.0 + float(lon_sec) / 3600.0
	lon_rads = lon_degs * pi / 180.0
	# changing azimuth of the source in degrees to radians
	az_deg = aa.__str__().split()[1].split('d')[0]
	az_min = aa.__str__().split()[2].split(''')[0]
	azsec = aa.__str__().split()[3].split(''')[0]
	az_degs = float(az_deg) + float(az_min) / 60.0 + float(azsec) / 3600.0
	az_rads = az_degs * pi / 180.0	
	# changing elevation or altitude of the source in degrees to radians
	al_deg = aa.__str__().split()[5].split('d')[0]
	al_min = aa.__str__().split()[6].split(''')[0]
	al_sec = aa.__str__().split()[7].split(''')[0]
	al_degs = float(al_deg) + float(al_min) / 60.0 + float(al_sec) / 3600.0
	al_rads = al_degs * pi / 180.0

	# all the values are returned in radians!
	return az_rads, al_rads, h, lat_rads, lon_rads


def check_args(ti):
	print('check_args')
	print(ti)
	'''
	Process all command line arguments.

	 [ if sys.argv[1:] is a valid set of command line arguments ->
		 return (ra_dec, lat_lon, dt) where ra_dec is a set of
		 celestial coordinates as a sidereal.ra_dec instance,
		 lat_lon is position as a sidereal.lat_lon instance, and
		 dt is a datetime.datetime instance
		else ->
		 sys.stderr +:= error message
		 stop execution ]
	'''
	#-- 1 --
	# [ if sys.argv[1:] has exactly four elements ->
	#	 raw_ra_dec, raw_lat, raw_lon, raw_dt := those elements
	# else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	arg_list = sys.argv[1:]
	if len(arg_list) != 5:
		usage ('Incorrect command line argument count.')
	else:
		raw_ra_dec, raw_lat, raw_lon, raw_dt, IONEX_TEC_file = arg_list
	raw_dt = str(ti)
	print('arg_list', arg_list)

	#-- 2 --
	# [ if raw_ra_dec is a valid set of equatorial coordinates ->
	#	 ra_dec := those coordinates as a sidereal.ra_dec instance
	# else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	ra_dec = check_ra_dec(raw_ra_dec)

	#-- 3 --
	# [ if raw_lat is a valid latitude ->
	#	 lat := that latitude in radians
	# else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	try:
		lat = sidereal.parse_lat(raw_lat)
	except SyntaxError, detail:
		usage('Invalid latitude: {detail}'.format(detail=detail))

	#-- 4 --
	# [ if raw_lon is a valid longitude ->
	#	 lon := that longitude in radians
	# else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	try:
		lon = sidereal.parse_lon(raw_lon)
	except SyntaxError, detail:
		usage('Invalid longitude: {detail}'.format(detail=detail))

	#-- 5 --
	# [ if raw_dt is a valid date-time string ->
	#	 dt := that date-time as a datetime.datetime instance
	# else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	try:
		dt = sidereal.parse_datetime(raw_dt)
	except SyntaxError, detail:
		usage('Invalid timestamp: {detail}'.format(detail=detail))

	#-- 6 --
	lat_lon = sidereal.lat_lon(lat, lon)
	return (ra_dec, lat_lon, dt)
#--- usage

def usage(*L):
	'''Print a usage message and stop.

	 [ L is a list of strings ->
		 sys.stderr +:= (usage message) + (elements of L,
						 concatenated)
		 stop execution ]
	'''
	print(>>sys.stderr, '*** Usage:')
	print(>>sys.stderr, '*** rdaa RA+dec lat lon datetime')
	print(>>sys.stderr, '*** Or:')
	print(>>sys.stderr, '*** rdaa RA-dec lat lon datetime')
	print(>>sys.stderr, '*** Error: {msg}'.format(msg=''.join(L)))
	raise SystemExit
#--- check_ra_dec

def check_ra_dec(raw_ra_dec):
	'''Check and convert a pair of equatorial coordinates.

	 [ raw_ra_dec is a string ->
		 if raw_ra_dec is a valid set of equatorial coordinates ->
			return those coordinates as a sidereal.ra_dec instance
		 else ->
			sys.stderr +:= error message
			stop execution ]
	'''
	#-- 1 --
	# [ if raw_ra_dec contains either a '+' or a '-' ->
	#	 m := a re.match instance describing the first matching
	#			character
	# else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	m = SIGN_PAT.search(raw_ra_dec)
	if m is None:
		usage("Equatorial coordinates must be separated by '+' or '-'.")

	#-- 2 --
	# [ raw_ra := raw_ra_dec up to the match described by m
	# sign := characters matched by m
	# raw_dec := raw_ra_dec past the match described by m ]
	raw_ra = raw_ra_dec[:m.start()]
	sign = m.group()
	raw_dec = raw_ra_dec[m.end():]

	#-- 3 --
	# [ if raw_ra is a valid hours expression ->
	#	 ra := raw_ra as radians
	# else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	try:
		ra_hours = sidereal.parse_hours(raw_ra)
		ra = sidereal.hours_to_radians(ra_hours)
	except SyntaxError, detail:
		usage("Right ascension '{raw_ra}' should have the form 'NNh[NNm[NN.NNNs]]'.".format(raw_ra=raw_ra))

	#-- 4 --
	# [ if raw_dec is a valid angle expression ->
	#	 abs_dec := that angle in radians
	#	 sys.stderr +:= error message
	#	 stop execution ]
	try:
		abs_dec = sidereal.parse_angle(raw_dec)
	except SyntaxError, detail:
		usage("Right ascension '{raw_ra}' should have the form 'NNd[NNm[NN.NNNs]]'.".format(raw_ra=raw_ra))
	#-- 5 --
	if sign == '-':
		dec = -abs_dec
	else:
		dec = abs_dec

	#-- 6 --
	return sidereal.ra_dec(ra, dec)
#================================================================
# Epilogue
#----------------------------------------------------------------

if __name__ == '__main__':
	main()
