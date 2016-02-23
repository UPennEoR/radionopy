#!/usr/bin/env python

#-------------------------------------------------------------------
# This function provides de geographic and topographic coordinates
# the Ionospheric piercing point (IPP)
# @version 1.0
# 
#
# Given the azimuth and zenith angle of the line of sight at the 
# location of the antenna, the offsets in geographic coordinates at 
# the intersection of the line of sight with the ionosphere are 
# calculated. Also, the altitude and azimuth corrdinates at the IPP
# are estimated. 
#
# The ionsphere is assumed to approximated by a thin shell at a
# uniform altitude.
#
# Input: 
#	lat_obs		latitude of the antenna (radians)
#	az_sou		Azimtuh of the source (radians)
#				from antenna location
#	ze_sou		Zenith of the source (radians)
#				from antenna location
#	alt_ion		height of the Ionospheric thin shell
#				(meters)
# Output: 
#	d_lat		offset latitude (radians)
#	d_lon		offset longitude (radians)
#	az_punc		Azimuth of the source (radians)
#				from IPP
#	zen_punc	Zenith of the source (radians)
#				from IPP
#-------------------------------------------------------------------

from scipy import *

def punc_ion_offset(lat_obs, az_sou, ze_sou, alt_ion):
	radius_earth = 6371000.0 # in meters

	# The 2-D sine rule gives the zenith angle at the
	# Ionospheric piercing point
	zen_punc = math.asin((radius_earth * math.sin(ze_sou)) / (radius_earth + alt_ion)) 

	# Use the sum of the internal angles of a triange to determine theta
	theta = ze_sou - zen_punc

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
