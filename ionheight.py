#!/usr/bin/env python

#------------------------------------------
# This script reads a IONEX file and retrieves
# the height of the Ionosphere
#------------------------------------------

import numpy

def calcionheight(filename): 
	# opening and reading the IONEX file into memory
	with open(filename, 'r') as read_file:
		linestring = read_file.read()
		LongList = linestring.split('\n')

	for file_data in LongList:
		if file_data.split()[-1] == 'DHGT'
			IonH = float(LongList[i].split()[0])

	return IonH



		
		




