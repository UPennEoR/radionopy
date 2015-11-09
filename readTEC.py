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
#	coordLat	latitude of the antenna (degrees)
#	coordLon	longitude of the antenna (degrees)
#	filename	IONEX file name
# Output: 
#	TEC		array containing TEC values
# 	TECvalues[LAT,LON] = [00,01,02,...,22,23,24]hrs
#------------------------------------------------------

from __future__ import print_function
import numpy as np

def readIonexTEC(filename):
	#==========================================================================
	# Reading and storing only the TEC values of 1 day
	# (13 maps) into a 3D array

	# Opening and reading the IONEX file into memory
	with open(filename, 'r') as read_file:
		linestring = read_file.read()
		LongList = linestring.split('\n')

	# creating a new array without the header and only
	# with the TEC maps
	add = 0 
	NewLongList = []
	for file_data in LongList:
		if file_data.split()[-2:] == ['RMS', 'MAP']:
			add = 0
		elif file_data.split()[-2:] == ['IN', 'FILE']:
			NumberOfMaps = float(file_data.split()[0])

		if add == 1:
			NewLongList.append(file_data)

		if file_data.split()[0] == 'END' and file_data.split()[2] == 'HEADER':
			add = 1

		if file_data.split()[-1] == 'DHGT':
			IonH = float(file_data.split()[0])
		elif file_data.split()[-1] == 'DLAT':
			startLat, endLat, stepLat = [float(data_item) for data_item in file_data[:2]]
		elif file_data.split()[-1] == 'DLON':
			startLon, endLon, stepLon = [float(data_item) for data_item in file_data[:2]]

	# Variables that indicate the number of points in Lat. and Lon.
	pointsLon = ((endLon - startLon) / stepLon) + 1
	pointsLat = ((endLat - startLat) / stepLat) + 1

	print(startLon, endLon, stepLon)
	print(startLat, endLat, stepLat)
	print(pointsLon, pointsLat)

	# What are the Lat/Lon coords?
	Longitude = np.linspace(startLon, endLon, num=pointsLon)
	Latitude = np.linspace(startLat, endLat, num=pointsLat)

	# 3D array that will contain TEC values only
	a = np.zeros((NumberOfMaps, pointsLat, pointsLon))

	# Selecting only the TEC values to store in the 3-D array
	counterMaps = 1
	for i in range(len(NewLongList)):
		# Pointing to first map (out of 13 maps)
		# then by changing 'counterMaps' the other
		# maps are selected
		if NewLongList[i].split()[0] == str(counterMaps) and NewLongList[i].split()[-4] == 'START':
			# pointing the starting Latitude
			# then by changing 'counterLat' we select
			# TEC data at other latitudes within
			# the selected map
			counterLat = 0
			newstartLat = float(str(startLat))
			for itemLat in range(int(pointsLat)):
				if NewLongList[i + 2 + counterLat].split()[0].split('-')[0] == str(newstartLat)\
				or '-' + NewLongList[i + 2 + counterLat].split()[0].split('-')[1] == str(newstartLat):
					# Adding to array 'a' a line of Latitude TEC data
					# we account for the TEC values at negative latitudes
					counterLon = 0
					for count_num in range(3, 8):
						list_index = i + count_num + counterLat
						for item in range(len(NewLongList[list_index].split())):
							a[counterMaps - 1, itemLat, counterLon] = NewLongList[list_index].split()[item]
							counterLon = counterLon + 1
				counterLat = counterLat + 6
				newstartLat = newstartLat + stepLat
			counterMaps = counterMaps + 1

	return {'TEC': np.array(a), 'lat': Latitude, 'lon': Longitude}
	#==========================================================================

# This was the original version as supplied in ionFR

def calcTEC(coordLat, coordLon, filename): 
	timeInt = 1.0 # hours
	totalmaps = 25

	#==========================================================================
	# Reading and storing only the TEC values of 1 day
	# (13 maps) into a 3D array

	# Opening and reading the IONEX file into memory
	with open(filename, 'r') as read_file:
		linestring = read_file.read()
		LongList = linestring.split('\n')

	# creating a new array without the header and only
	# with the TEC maps
	add = 0 
	NewLongList = []
	for file_data in LongList:
		if file_data.split()[-2:] == ['RMS', 'MAP']:
			add = 0
		elif file_data.split()[-2:] == ['IN', 'FILE']:
			NumberOfMaps = float(file_data.split()[0])

		if add == 1:
			NewLongList.append(file_data)

		if file_data.split()[0] == 'END' and file_data.split()[2] == 'HEADER':
			add = 1

		if file_data.split()[-1] == 'DHGT':
			IonH = float(file_data.split()[0])
		elif file_data.split()[-1] == 'DLAT':
			startLat, endLat, stepLat = [float(data_item) for data_item in file_data[:2]]
		elif file_data.split()[-1] == 'DLON':
			startLon, endLon, stepLon = [float(data_item) for data_item in file_data[:2]]

	# Variables that indicate the number of points in Lat. and Lon.
	pointsLon = ((endLon - startLon) / stepLon) + 1
	pointsLat = ((endLat - startLat) / stepLat) + 1

	# 3D array that will contain TEC values only
	a = np.zeros((NumberOfMaps, pointsLat, pointsLon))

	# Selecting only the TEC values to store in the 3-D array
	counterMaps = 1
	for i in range(len(NewLongList)):
		# Pointing to first map (out of 13 maps)
		# then by changing 'counterMaps' the other
		# maps are selected
		if NewLongList[i].split()[0] == str(counterMaps) and NewLongList[i].split()[-4] == 'START':
			# pointing the starting Latitude
			# then by changing 'counterLat' we select
			# TEC data at other latitudes within
			# the selected map
			counterLat = 0
			newstartLat = float(str(startLat))
			for itemLat in range(int(pointsLat)):
				if NewLongList[i + 2 + counterLat].split()[0].split('-')[0] == str(newstartLat)\
				or '-' + NewLongList[i + 2 + counterLat].split()[0].split('-')[1] == str(newstartLat):
					# Adding to array 'a' a line of Latitude TEC data
					# we account for the TEC values at negative latitudes
					counterLon = 0
					for count_num in range(3, 8):
						list_index = i + count_num + counterLat
						for item in range(len(NewLongList[list_index].split())):
							a[counterMaps - 1, itemLat, counterLon] = NewLongList[list_index].split()[item]
							counterLon = counterLon + 1
				counterLat = counterLat + 6
				newstartLat = newstartLat + stepLat
			counterMaps = counterMaps + 1
	#==========================================================================


	#==========================================================================================
	# producing interpolated TEC maps, and consequently a new array that will 
	# contain 25 TEC maps in total. The interpolation method used is the second
	# one indicated in the IONEX manual

	# creating a new array that will contain 25 maps in total 
	newa = np.zeros((totalmaps, pointsLat, pointsLon))
	inc = 0
	for item in range(int(NumberOfMaps)):
		newa[inc, :, :] = a[item, :, :]
		inc = inc + 2

	# performing the interpolation to create 12 addional maps 
	# from the 13 TEC maps available
	while int(timeInt) <= (totalmaps - 2):
		for lat in range(int(pointsLat)):
			for lon in range(int(pointsLon)):
				# interpolation type 2:
				# newa[int(timeInt),lat,lon] = 0.5*newa[int(timeInt)-1,lat,lon] + 0.5*newa[int(timeInt)+1,lat,lon]
				# interpolation type 3 ( 3 or 4 columns to the right and left of the odd maps have values of zero
				# Correct for this):
				if (lon >= 4) and (lon <= (pointsLon - 4)):
					newa[int(timeInt), lat, lon] = 0.5 * newa[int(timeInt) - 1, lat, lon + 3] + 0.5 * newa[int(timeInt) + 1, lat, lon -3 ] 
		timeInt = timeInt + 2.0
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
	for lon in range(int(pointsLon)):
		if (coordLon > (startLon + (n + 1) * stepLon)  and coordLon < (startLon + (n + 2) * stepLon)):
			lowerIndexLon =  n + 1
			higherIndexLon = n + 2
		n = n + 1
	for lat in range(int(pointsLat)):
		if (coordLat < (startLat + (m + 1) * stepLat)  and coordLat > (startLat + (m + 2) * stepLat)):
			lowerIndexLat =  m + 1
			higherIndexLat = m + 2
		m = m + 1

	# Using the 4-point formula indicated in the IONEX manual
	# The TEC value at the coordinates you desire for every 
	# hour are estimated 
	diffLon = coordLon - (startLon + lowerIndexLon * stepLon)
	p = diffLon / stepLon
	diffLat = coordLat - (startLat + lowerIndexLat * stepLat)
	q = diffLat / stepLat
	TECvalues = []
	for m in range(totalmaps):
		TECvalues.append((1.0 - p) * (1.0 - q) * newa[m, lowerIndexLat, lowerIndexLon]\
							+ p * (1.0 - q) * newa[m, lowerIndexLat, higherIndexLon]\
							+ q * (1.0 - p) * newa[m, higherIndexLat, lowerIndexLon]\
							+ p * q * newa[m, higherIndexLat, higherIndexLon])
	#=========================================================================

	return {'TECvalues': np.array(TECvalues), 'a': np.array(a), 'newa': np.array(newa)}
