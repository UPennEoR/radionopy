#!/usr/bin/env python

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
