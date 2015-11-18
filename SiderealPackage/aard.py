#!/usr/bin/env python
from __future__ import print_function
import sys
import re
import sidereal as sdr
#================================================================
# Manifest consants
#----------------------------------------------------------------

SIGN_PAT = re.compile(r'[\-+]')
#-----main

def main():
	'''
	Main program for aard.
	'''

	#-- 1 --
	# [ if sys.argv contains a valid set of command line
	#  arguments ->
	#	 alt_az := the azimuth and altitude as
	#				a sdr.alt_az instance
	#	 lat_lon := the observer's location as a
	#				 sdr.lat_lon instance
	#	 dt := the observer's date and time as a
	#			 datetime.datetime instance
	#  else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	alt_az, lat_lon, dt = check_args()

	#-- 2 --
	# [ if dt has no time zone information ->
	#	 utc := dt
	#  else ->
	#	 utc := the UTC equivalent to dt ]
	if dt.tzinfo is None or dt.utc_offset() is None:
		utc = dt
	else:
		utc = dt - dt.utc_offset()

	#-- 3 --
	# [ sys.stdout +:= local sidereal time for dt and lat_lon ]
	gst = sdr.SiderealTime.from_datetime(utc)
	lst = gst.lst(lat_lon.lon)
	print('Horizon coordinates:', alt_az)
	print('Observer\'s location:', lat_lon)
	print('Observer\'s time:', dt)
	print('Local sidereal time is', lst)

	#-- 4 --
	# [ ra_dec := equatorial coordinates of self for local
	#	  sidereal time (lst) and location (lat_lon) ]
	ra_dec = alt_az.ra_dec(lst, lat_lon)

	#-- 5 --
	print('Equatorial coordinates:', ra_dec)
#--- check_args

def check_args():
	'''
	Process all command line arguments.

	 [ if sys.argv[1:] is a valid set of command line arguments ->
		 return (alt_az, lat_lon, dt) where alt_az is a set of
		 horizon coordinates as a sdr.alt_az instance,
		 lat_lon is position as a sdr.lat_lon instance, and
		 dt is a datetime.datetime instance
		else ->
		 sys.stderr +:= error message
		 stop execution ]
	'''

	#-- 1 --
	# [ if sys.argv[1:] has exactly four elements ->
	#	 raw_alt_az, raw_lat, raw_lon, raw_dt := those elements
	#  else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	arg_list = sys.argv[1:]
	if len(arg_list) != 4:
		usage ('Incorrect command line argument count.')
	else:
		raw_alt_az, raw_lat, raw_lon, raw_dt = arg_list

	#-- 2 --
	# [ if raw_alt_az is a valid set of horizon coordinates ->
	#	 alt_az := those coordinates as a sdr.alt_az instance
	alt_az = check_alt_az(raw_alt_az)

	#-- 3 --
	# [ if raw_lat is a valid latitude ->
	#	 lat := that latitude in radians
	#  else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	try:
		lat = sdr.parse_lat(raw_lat)
	except SyntaxError, detail:
		usage('Invalid latitude: {detail}'.format(detail=detail))

	#-- 4 --
	# [ if raw_lon is a valid longitude ->
	#	 lon := that longitude in radians
	#  else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	try:
		lon = sdr.parse_lon(raw_lon)
	except SyntaxError, detail:
		usage('Invalid longitude: {detail}'.format(detail=detail))

	#-- 5 --
	# [ if raw_dt is a valid date-time string ->
	#	 dt := that date-time as a datetime.datetime instance
	#  else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	try:
		dt = sdr.parse_datetime(raw_dt)
	except SyntaxError, detail:
		usage('Invalid timestamp: {detail}'.format(detail=detail))

	#-- 6 --
	lat_lon = sdr.lat_lon(lat, lon)
	return (alt_az, lat_lon, dt)
# - - -  u s a g e

def usage(*L):
	'''Print a usage message and stop.

	 [ L is a list of strings ->
		 sys.stderr +:= (usage message) + (elements of L,
						  concatenated)
		 stop execution ]
	'''
	print(>>sys.stderr, '*** Usage:')
	print(>>sys.stderr, '***  aard az+alt lat lon datetime')
	print(>>sys.stderr, '*** Error: {msg}'.format(msg=''.join(L)))
	raise SystemExit
#--- check_alt_az

def check_alt_az(raw_alt_az):
	'''
	Check and convert a pair of horizon coordinates.

	 [ raw_alt_az is a string ->
		 if raw_alt_az is a valid set of horizon coordinates ->
			return those coordinates as a sdr.alt_az instance
		 else ->
			sys.stderr +:= error message
			stop execution ]
	'''
	#-- 1 --
	# [ if raw_alt_az contains either a '+' or a '-' ->
	#	 m := a re.match instance describing the first matching
	#			character
	#  else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	m = SIGN_PAT.search(raw_alt_az)
	if m is None:
		usage("Equatorial coordinates must be separated by '+' or '-'.")

	#-- 2 --
	# [ raw_az := raw_alt_az up to the match described by m
	#  sign := characters matched by m
	#  raw_alt := raw_alt_az past the match described by m ]
	raw_az = raw_alt_az[:m.start()]
	sign = m.group()
	raw_alt = raw_alt_az[m.end():]

	#-- 3 --
	# [ if raw_az is a valid angle ->
	#	 az := that angle as radians
	#  else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	try:
		az = sdr.parse_angle(raw_az)
	except SyntaxError, detail:
		usage("Azimuth '{az}' should have the form 'NNNd[NNm[NN.NNNs]]'.".format(az=raw_az))

	#-- 4 --
	# [ if raw_alt is a valid angle ->
	#	 alt := that angle as radians
	#  else ->
	#	 sys.stderr +:= error message
	#	 stop execution ]
	try:
		abs_alt = sdr.parse_angle(raw_alt)
	except SyntaxError, detail:
		usage("Altitude '{alt}' should have the form 'NNd[NNm[NN.NNNs]]'.".format(alt=raw_alt)

	#-- 5 --
	if sign == '-':
		alt = -abs_alt
	else:
		alt = abs_alt

	#-- 6 --
	return sdr.alt_az(alt, az)
#================================================================
# Epilogue
#----------------------------------------------------------------

if __name__ == '__main__':
	main()
