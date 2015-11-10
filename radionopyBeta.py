import ftplib
import datetime
import os
import numpy as np
from astropy import units as u
from astropy import constants as c


# From ionFRM.py
# Defining some variables for further use
TECU = pow(10,16)
TEC2m2 = 0.1*TECU
EarthRadius = 6371000.0 # in meters
# This differs a bit from what ionFRM uses, by about 7 km
RadiusEarth = c.R_earth.value 
Tesla2Gauss = np.power(10,4)

# ------------------------------------------------------------------------------
# Generate the necessary IONEX filename for a given day, and fetch it via ftp
# ------------------------------------------------------------------------------
# From the little script from ionFR of the same name
def IONEXFileNeeded(year,month,day):
    dayofyear = datetime.datetime.strptime(''+str(year)+' '+str(month)+' '+str(day)+'', '%Y %m %d').timetuple().tm_yday

    if dayofyear < 10:
        dayofyear = '00'+str(dayofyear)
    if dayofyear < 100 and dayofyear >= 10:
        dayofyear = '0'+str(dayofyear)

    # Outputing the name of the IONEX file you require
    file =  'CODG'+str(dayofyear)+'0.'+str(list(str(year))[2])+''+str(list(str(year))[3])+'I'
    return file

server = 'ftp.unibe.ch'

def getIONEXfile(year,month,day):

    ftp_dir = 'aiub/CODE/'+year
    ionexfile = IONEXFileNeeded(year,month,day)
    ionexfileZ = ionexfile+'.Z'

    print 'Retrieving '+ionexfileZ+' for '+day+' '+month+' '+year
    ftp = ftplib.FTP(server, 'anonymous', 'jaguirre@sas.upenn.edu')
    ftp.cwd(ftp_dir)
    ftp.retrbinary('RETR '+ionexfileZ, open(ionexfileZ, 'wb').write)
    ftp.quit()

    return True

#    os.system('gunzip '+ionexfileZ)

#------------------------------------------------------
# Extract TEC values from an IONEX file
# given a specific geographic coordinate.
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

def readIonexTEC(filename):

#==========================================================================
	# Reading and storing only the TEC values of 1 day
	# (13 maps) into a 3D array

	# Opening and reading the IONEX file into memory
	linestring = open(filename, 'r').read()
	LongList = linestring.split('\n')

	# creating a new array without the header and only
	# with the TEC maps
	add = 0 
	NewLongList = []
	for i in range(len(LongList)-1):
                if LongList[i].split()[-1] == 'DHGT':
			IonH = float(LongList[i].split()[0])
		if LongList[i].split()[-1] == 'MAP':
			if LongList[i].split()[-2] == 'RMS':
				add = 0
		if add == 1:	
			NewLongList.append(LongList[i])
		if LongList[i].split()[-1] == 'FILE':
			if LongList[i].split()[-2] == 'IN':
				NumberOfMaps = float(LongList[i].split()[0])
		if LongList[i].split()[-1] == 'DHGT':
			IonH = float(LongList[i].split()[0])
		if LongList[i].split()[-1] == 'DLAT':
			startLat = float(LongList[i].split()[0])
			endLat = float(LongList[i].split()[1])
			stepLat = float(LongList[i].split()[2])
		if LongList[i].split()[-1] == 'DLON':
			startLon = float(LongList[i].split()[0])
			endLon = float(LongList[i].split()[1])
			stepLon = float(LongList[i].split()[2])
		if LongList[i].split()[0] == 'END':
			if LongList[i].split()[2] == 'HEADER':
				add = 1	

	# Variables that indicate the number of points in Lat. and Lon.
	pointsLon = ((endLon - startLon)/stepLon) + 1
	pointsLat = ((endLat - startLat)/stepLat) + 1

        print startLon,endLon,stepLon
        print startLat,endLat,stepLat
        print pointsLon,pointsLat
        
        # What are the Lat/Lon coords?
        Longitude = np.linspace(startLon,endLon,num=pointsLon)
        Latitude = np.linspace(startLat,endLat,num=pointsLat)
        
	# 3D array that will contain TEC values only
	a = np.zeros((NumberOfMaps,pointsLat,pointsLon))

	# Selecting only the TEC values to store in the 3-D array
	counterMaps = 1
	for i in range(len(NewLongList)):
		# Pointing to first map (out of 13 maps)
		# then by changing 'counterMaps' the other
		# maps are selected
		if NewLongList[i].split()[0] == ''+str(counterMaps)+'':
			if NewLongList[i].split()[-4] == 'START':
				# pointing the starting Latitude
				# then by changing 'counterLat' we select
				# TEC data at other latitudes within
				# the selected map
				counterLat = 0
				newstartLat = float(str(startLat))
				for itemLat in range(int(pointsLat)):
					if NewLongList[i+2+counterLat].split()[0].split('-')[0] == ''+str(newstartLat)+'':
						# adding to array 'a' a line of Latitude TEC data
						counterLon = 0
						for item in range(len(NewLongList[i+3+counterLat].split())):
							a[counterMaps-1,itemLat,counterLon] = NewLongList[i+3+counterLat].split()[item]
							counterLon = counterLon + 1
						for item in range(len(NewLongList[i+4+counterLat].split())):
							a[counterMaps-1,itemLat,counterLon] = NewLongList[i+4+counterLat].split()[item]
							counterLon = counterLon + 1
						for item in range(len(NewLongList[i+5+counterLat].split())):
							a[counterMaps-1,itemLat,counterLon] = NewLongList[i+5+counterLat].split()[item]
							counterLon = counterLon + 1
						for item in range(len(NewLongList[i+6+counterLat].split())):
							a[counterMaps-1,itemLat,counterLon] = NewLongList[i+6+counterLat].split()[item]
							counterLon = counterLon + 1
						for item in range(len(NewLongList[i+7+counterLat].split())):
							a[counterMaps-1,itemLat,counterLon] = NewLongList[i+7+counterLat].split()[item]
							counterLon = counterLon + 1				
					if '-'+NewLongList[i+2+counterLat].split()[0].split('-')[1] == ''+str(newstartLat)+'':
						# Adding to array 'a' a line of Latitude TEC data
						# Same chunk as above but in this case we account for
						# the TEC values at negative latitudes
						counterLon = 0
						for item in range(len(NewLongList[i+3+counterLat].split())):
							a[counterMaps-1,itemLat,counterLon] = NewLongList[i+3+counterLat].split()[item]
							counterLon = counterLon + 1
						for item in range(len(NewLongList[i+4+counterLat].split())):
							a[counterMaps-1,itemLat,counterLon] = NewLongList[i+4+counterLat].split()[item]
							counterLon = counterLon + 1
						for item in range(len(NewLongList[i+5+counterLat].split())):
							a[counterMaps-1,itemLat,counterLon] = NewLongList[i+5+counterLat].split()[item]
							counterLon = counterLon + 1
						for item in range(len(NewLongList[i+6+counterLat].split())):
							a[counterMaps-1,itemLat,counterLon] = NewLongList[i+6+counterLat].split()[item]
							counterLon = counterLon + 1
						for item in range(len(NewLongList[i+7+counterLat].split())):
							a[counterMaps-1,itemLat,counterLon] = NewLongList[i+7+counterLat].split()[item]
							counterLon = counterLon + 1
					counterLat = counterLat + 6
					newstartLat = newstartLat + stepLat
				counterMaps = counterMaps + 1

        return {'TEC':np.array(a),'lat':Latitude,'lon':Longitude,'AltIon':IonH*1000.0}

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
#	LatObs		latitude of the antenna (radians)
#	AzSou		Azimtuh of the source (radians)
#			from antenna location
#	ZeSou		Zenith of the source (radians)
#			from antenna location
#	AltIon		height of the Ionospheric thin shell
#			(meters)
# Output: 
#	dLat		offset latitude (radians)
#	dLon		offset longitude (radians)
#	AzPunc		Azimuth of the source (radians)
#			from IPP
#	ZenPunc		Zenith of the source (radians)
#			from IPP
#-------------------------------------------------------------------

#from scipy import *

def PuncIonOffset(LatObs,AzSou,ZeSou,AltIon):

	#RadiusEarth = 6371000.0 # in meters

	# The 2-D sine rule gives the zenith angle at the
	# Ionospheric piercing point
	ZenPunc = np.arcsin((RadiusEarth*np.sin(ZeSou))/(RadiusEarth + AltIon)) 

	# Use the sum of the internal angles of a triange to determine theta
	theta = ZeSou - ZenPunc

	# The cosine rule for spherical triangles gives us the latitude
	# at the IPP
	lation = np.arcsin(np.sin(LatObs)*np.cos(theta) + np.cos(LatObs)*np.sin(theta)*np.cos(AzSou)) 
	dLat = lation - LatObs # latitude difference

	# Longitude difference using the 3-D sine rule (or for spherical triangles)
	dLon = np.arcsin(np.sin(AzSou)*np.sin(theta)/np.cos(lation))

	# Azimuth at the IPP using the 3-D sine rule
	sazion = np.sin(AzSou)*np.cos(LatObs)/np.cos(lation)
	AzPunc = np.arcsin(sazion)

	return dLat,dLon,AzPunc,ZenPunc

#------------------------------------------
# This script reads a IONEX file and retrieves
# the height of the Ionosphere
#------------------------------------------

#import numpy

def calcionheight(filename): 

	# opening and reading the IONEX file into memory
	linestring = open(filename, 'r').read()
	LongList = linestring.split('\n')
	################################################

	for i in range(len(LongList)-1):
		if LongList[i].split()[-1] == 'DHGT':
			IonH = float(LongList[i].split()[0])

	return IonH



		
		






		
		





