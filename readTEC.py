#!/usr/bin/env python

#------------------------------------------------------
# Extract TEC values from an IONEX file
# given a specific geographic coordinate.
# @version 1.0
# @author carlos
#
# The TEC measurements provided in the IONEX
# files are vertical TEC values, i.e. TEC values
# at the zenith. The TEC values you eventually desire
# to use for the computation of the IFR have to be
# divided by the cos(ZenithSource -> the direction
# along the line of sight of the source of interest).
#
# 25 TEC maps from the 13 initially provided are
# created. The interpolation method used is the third
# one indicated in the IONEX manual. A grid interpolation
# is also used to find out the 'exact' TEC value
# at the coordinates you require.
#
# Input: 
#	coord_lat	latitude of the antenna (degrees)
#	coord_lon	longitude of the antenna (degrees)
#	filename	IONEX file name
# Output: 
#	TEC		array containing TEC values
# 	TECvalues[LAT,LON] = [00,01,02,...,22,23,24]hrs
#------------------------------------------------------

from __future__ import print_function
import numpy as np

def read_IONEX_TEC(filename):
	#==========================================================================
	# Reading and storing only the TEC values of 1 day
	# (13 maps) into a 3D array

	# Opening and reading the IONEX file into memory
	with open(filename, 'r') as read_file:
		linestring = read_file.read()
		IONEX_list = linestring.split('\n')

	# creating a new array without the header and only
	# with the TEC maps
	add = 0 
	new_IONEX_list = []
	for file_data in IONEX_list:
		if file_data.split()[-2:] == ['RMS', 'MAP']:
			add = 0
		elif file_data.split()[-2:] == ['IN', 'FILE']:
			number_of_maps = float(file_data.split()[0])

		if add == 1:
			new_IONEX_list.append(file_data)

		if file_data.split()[0] == 'END' and file_data.split()[2] == 'HEADER':
			add = 1

		if file_data.split()[-1] == 'DHGT':
			ion_h = float(file_data.split()[0])
		elif file_data.split()[-1] == 'DLAT':
			start_lat, end_lat, step_lat = [float(data_item) for data_item in file_data[:2]]
		elif file_data.split()[-1] == 'DLON':
			start_lon, end_lon, step_lon = [float(data_item) for data_item in file_data[:2]]

	# Variables that indicate the number of points in Lat. and Lon.
	points_lon = ((end_lon - start_lon) / step_lon) + 1
	points_lat = ((end_lat - start_lat) / step_lat) + 1

	print(start_lon, end_lon, step_lon)
	print(start_lat, end_lat, step_lat)
	print(points_lon, points_lat)

	# What are the Lat/Lon coords?
	longitude = np.linspace(start_lon, end_lon, num=points_lon)
	latitude = np.linspace(start_lat, end_lat, num=points_lat)

	# 3D array that will contain TEC values only
	a = np.zeros((number_of_maps, points_lat, points_lon))

	# Selecting only the TEC values to store in the 3-D array
	counter_maps = 1
	for i in range(len(new_IONEX_list)):
		# Pointing to first map (out of 13 maps)
		# then by changing 'counter_maps' the other
		# maps are selected
		if new_IONEX_list[i].split()[0] == str(counter_maps) and new_IONEX_list[i].split()[-4] == 'START':
			# pointing the starting latitude
			# then by changing 'counter_lat' we select
			# TEC data at other latitudes within
			# the selected map
			counter_lat = 0
			new_start_lat = float(str(start_lat))
			for item_lat in range(int(points_lat)):
				if new_IONEX_list[i + 2 + counter_lat].split()[0].split('-')[0] == str(new_start_lat)\
				or '-' + new_IONEX_list[i + 2 + counter_lat].split()[0].split('-')[1] == str(new_start_lat):
					# Adding to array 'a' a line of latitude TEC data
					# we account for the TEC values at negative latitudes
					counter_lon = 0
					for count_num in range(3, 8):
						list_index = i + count_num + counter_lat
						for new_IONEX_item in new_IONEX_list[list_index].split():
							a[counter_maps - 1, item_lat, counter_lon] = new_IONEX_item
							counter_lon = counter_lon + 1
				counter_lat = counter_lat + 6
				new_start_lat = new_start_lat + step_lat
			counter_maps = counter_maps + 1

	return {'TEC': np.array(a), 'lat': latitude, 'lon': longitude}
	#==========================================================================

# This was the original version as supplied in ionFR

def calc_TEC(coord_lat, coord_lon, filename): 
	time_int = 1.0 # hours
	total_maps = 25

	#==========================================================================
	# Reading and storing only the TEC values of 1 day
	# (13 maps) into a 3D array

	# Opening and reading the IONEX file into memory
	with open(filename, 'r') as read_file:
		linestring = read_file.read()
		IONEX_list = linestring.split('\n')

	# creating a new array without the header and only
	# with the TEC maps
	add = 0 
	new_IONEX_list = []
	for file_data in IONEX_list:
		if file_data.split()[-2:] == ['RMS', 'MAP']:
			add = 0
		elif file_data.split()[-2:] == ['IN', 'FILE']:
			number_of_maps = float(file_data.split()[0])

		if add == 1:
			new_IONEX_list.append(file_data)

		if file_data.split()[0] == 'END' and file_data.split()[2] == 'HEADER':
			add = 1

		if file_data.split()[-1] == 'DHGT':
			ion_h = float(file_data.split()[0])
		elif file_data.split()[-1] == 'DLAT':
			start_lat, end_lat, step_lat = [float(data_item) for data_item in file_data[:2]]
		elif file_data.split()[-1] == 'DLON':
			start_lon, end_lon, step_lon = [float(data_item) for data_item in file_data[:2]]

	# Variables that indicate the number of points in Lat. and Lon.
	points_lon = ((end_lon - start_lon) / step_lon) + 1
	points_lat = ((end_lat - start_lat) / step_lat) + 1

	# 3D array that will contain TEC values only
	a = np.zeros((number_of_maps, points_lat, points_lon))

	# Selecting only the TEC values to store in the 3-D array
	counter_maps = 1
	for i in range(len(new_IONEX_list)):
		# Pointing to first map (out of 13 maps)
		# then by changing 'counter_maps' the other
		# maps are selected
		if new_IONEX_list[i].split()[0] == str(counter_maps) and new_IONEX_list[i].split()[-4] == 'START':
			# pointing the starting latitude
			# then by changing 'counter_lat' we select
			# TEC data at other latitudes within
			# the selected map
			counter_lat = 0
			new_start_lat = float(str(start_lat))
			for item_lat in range(int(points_lat)):
				if new_IONEX_list[i + 2 + counter_lat].split()[0].split('-')[0] == str(new_start_lat)\
				or '-' + new_IONEX_list[i + 2 + counter_lat].split()[0].split('-')[1] == str(new_start_lat):
					# Adding to array 'a' a line of latitude TEC data
					# we account for the TEC values at negative latitudes
					counter_lon = 0
					for count_num in range(3, 8):
						list_index = i + count_num + counter_lat
						for new_IONEX_item in new_IONEX_list[list_index].split():
							a[counter_maps - 1, item_lat, counter_lon] = new_IONEX_item
							counter_lon = counter_lon + 1
				counter_lat = counter_lat + 6
				new_start_lat = new_start_lat + step_lat
			counter_maps = counter_maps + 1
	#==========================================================================


	#==========================================================================================
	# producing interpolated TEC maps, and consequently a new array that will 
	# contain 25 TEC maps in total. The interpolation method used is the second
	# one indicated in the IONEX manual

	# creating a new array that will contain 25 maps in total 
	newa = np.zeros((total_maps, points_lat, points_lon))
	inc = 0
	for item in range(int(number_of_maps)):
		newa[inc, :, :] = a[item, :, :]
		inc = inc + 2

	# performing the interpolation to create 12 addional maps 
	# from the 13 TEC maps available
	while int(time_int) <= (total_maps - 2):
		for lat in range(int(points_lat)):
			for lon in range(int(points_lon)):
				# interpolation type 2:
				# newa[int(time_int),lat,lon] = 0.5*newa[int(time_int)-1,lat,lon] + 0.5*newa[int(time_int)+1,lat,lon]
				# interpolation type 3 ( 3 or 4 columns to the right and left of the odd maps have values of zero
				# Correct for this):
				if (lon >= 4) and (lon <= (points_lon - 4)):
					newa[int(time_int), lat, lon] = 0.5 * newa[int(time_int) - 1, lat, lon + 3] + 0.5 * newa[int(time_int) + 1, lat, lon -3 ] 
		time_int = time_int + 2.0
	#==========================================================================================


	#=========================================================================
	# Finding out the TEC value for the coordinates given
	# at every hour

	# Locating the 4 points in the IONEX grid map which surround
	# the coordinate you want to calculate the TEC value from  
	indexLat = 0
	indexLon = 0
	n = 0
	m = 0
	for lon in range(int(points_lon)):
		if (coord_lon > (start_lon + (n + 1) * step_lon)  and coord_lon < (start_lon + (n + 2) * step_lon)):
			lower_index_lon =  n + 1
			higher_index_lon = n + 2
		n = n + 1
	for lat in range(int(points_lat)):
		if (coord_lat < (start_lat + (m + 1) * step_lat)  and coord_lat > (start_lat + (m + 2) * step_lat)):
			lower_index_lat =  m + 1
			higher_index_lat = m + 2
		m = m + 1

	# Using the 4-point formula indicated in the IONEX manual
	# The TEC value at the coordinates you desire for every 
	# hour are estimated 
	diff_lon = coord_lon - (start_lon + lower_index_lon * step_lon)
	p = diff_lon / step_lon
	diff_lat = coord_lat - (start_lat + lower_index_lat * step_lat)
	q = diff_lat / step_lat
	TECvalues = []
	for m in range(total_maps):
		TECvalues.append((1.0 - p) * (1.0 - q) * newa[m, lower_index_lat, lower_index_lon]\
							+ p * (1.0 - q) * newa[m, lower_index_lat, higher_index_lon]\
							+ q * (1.0 - p) * newa[m, higher_index_lat, lower_index_lon]\
							+ p * q * newa[m, higher_index_lat, higher_index_lon])
	#=========================================================================

	return {'TECvalues': np.array(TECvalues), 'a': np.array(a), 'newa': np.array(newa)}
