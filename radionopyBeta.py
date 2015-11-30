import ftplib
import datetime
import os
import numpy as np
from astropy import units as u
from astropy import constants as c
import sys

# From ionFRM.py
# Defining some variables for further use
TECU = pow(10,16)
TEC2m2 = 0.1*TECU
EarthRadius = 6371000.0 # in meters
# This differs a bit from what ionFRM uses, by about 7 km
RadiusEarth = c.R_earth.value 
Tesla2Gauss = np.power(10,4)

# The IONEX FTP server
server = 'ftp.unibe.ch'

# ------------------------------------------------------------------------------
# Generate the necessary IONEX filename for a given day, and fetch it via ftp
# ------------------------------------------------------------------------------
# From the little script from ionFR of the same name
def IONEXFileNeeded(year,month,day):
    dayofyear = datetime.datetime.strptime(''+str(year)+' '+str(month)+' '+str(day)+'', '%Y %m %d').timetuple().tm_yday

#    if dayofyear < 10:
#        dayofyear = '00'+str(dayofyear)
#    if dayofyear < 100 and dayofyear >= 10:
#        dayofyear = '0'+str(dayofyear)
    dayofyear=str(dayofyear).zfill(4)

    # Outputing the name of the IONEX file you require
    file =  'CODG'+str(dayofyear)+'0.'+str(list(str(year))[2])+''+str(list(str(year))[3])+'I'
    return file

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

# I can't believe how complicated this guy made this
def readIONEX(filename):
    # Read in the file
    FullFile = open(filename, 'r').read()
    header, data = FullFile.split('END OF HEADER')
    header = header.split('\n')
    data = data.split('\n')
    # Flip through the header and pull out useful information
    for h in header:
        if ('# OF MAPS IN FILE' in h):
            Nmaps = int(h.split()[0])
        if ('HGT1 / HGT2 / DHGT' in h):
            ionosphere_height = float(h.split()[0])
        # Get the size and coordinates of the map in latitude
        if ('LAT1 / LAT2 / DLAT' in h):
            lat1,lat2,dlat = [float(i) for i in h.split()[0:3]]
            lat = np.arange(lat1,lat2+dlat,dlat)
            Nlat = len(lat)
        # Get the size and coordinates of the map in longitude
        if ('LON1 / LON2 / DLON' in h):
            lon1,lon2,dlon = [float(i) for i in h.split()[0:3]]
            lon = np.arange(lon1,lon2+dlon,dlon)
            Nlon = len(lon)

    # Allocate the map arrays
    TEC = np.zeros([Nmaps,Nlat,Nlon])
    TEC_RMS = np.zeros([Nmaps,Nlat,Nlon])

    tec_start_indx = []
    tec_end_indx = []
    rms_start_indx = []
    rms_end_indx = []
    
    # Find the boundaries
    for i,d in enumerate(data):
        if ('START OF TEC MAP' in d):
            tec_start_indx.append(i)
        if ('END OF TEC MAP' in d):
            tec_end_indx.append(i)
        if ('START OF RMS MAP' in d):
            rms_start_indx.append(i)
        if ('END OF RMS MAP' in d):
            rms_end_indx.append(i)

    if (len(tec_start_indx) != Nmaps or len(tec_end_indx) != Nmaps):
        print 'Error: not the right number of starts and ends'

    for map in np.arange(Nmaps):
        # Build the TEC maps
        lines = data[tec_start_indx[map]:tec_end_indx[map]]
        vals = []
        for l in lines:
            if (('START OF TEC MAP' in l) or ('LAT/LON1/LON2/DLON/H' in l) or ('EPOCH OF CURRENT MAP' in l)):
                continue
            else:
                for v in l.split():
                    vals.append(float(v))
        TEC[map,:,:] = np.array(vals).reshape([Nlat,Nlon])

        # Build the RMS maps
        lines = data[rms_start_indx[map]:rms_end_indx[map]]
        vals = []
        for l in lines:
            if (('START OF RMS MAP' in l) or ('LAT/LON1/LON2/DLON/H' in l) or ('EPOCH OF CURRENT MAP' in l)):
                continue
            else:
                for v in l.split():
                    vals.append(float(v))
        TEC_RMS[map,:,:] = np.array(vals).reshape([Nlat,Nlon])
       
    return {'header':header,'data_in_file':data,'Nmaps':Nmaps,'IonHeight':ionosphere_height,'lat1':lat1,'lat2':lat2,'dlat':dlat,'lat':lat,'Nlat':Nlat,'lon1':lon1,'lon2':lon2,'dlon':dlon,'lon':lon,'Nlon':Nlon,'tec_i':tec_start_indx,'tec_f':tec_end_indx,'rms_i':rms_start_indx,'rms_f':rms_end_indx,'TEC':TEC,'TEC_RMS':TEC_RMS}
    

    
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
        # This is Equation 5
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

#-----------------------------------------------------------    
path='/Users/jaguirre/PyModules/'
# Maybe not necessary ...
#sys.path.append(""+str(path)+"ionFR/SiderealPackage")
#sys.path.append(""+str(path)+"ionFR/PunctureIonosphereCoord")
#sys.path.append(""+str(path)+"ionFR/IONEX")
#-----------------------------------------------------------

# wrap up the IGRF call
def B_IGRF(year,month,day,AltIon,LonO,offLon,LatO,offLat,AzPunct,ZenPunct):
    # Calculation of the total magnetic field along the line of sight at the IPP

    # Open the input file to geomag70.exe
    IGRFpath = ''+str(path)+'radionopy/IGRF/geomag70_linux/'
    f = open(IGRFpath+'input.txt', 'w')

    datestr = ''+str(year)+','+str(month)+','+str(day)
    IonRadiusstr = str((EarthRadius+AltIon)/1000.0)
    # Why does he bother with the original E/W N/S Lat/Lon entry?  It just gets parsed out everywhere ...
    # Ugh.  The minus sign is applied to (Lon + offLon), which means the sense of the offset is different, depending on the hemisphere
    # Are you FUCKING with me?
    Latstr = str(np.degrees(LatO + offLat))
    Lonstr = str(np.degrees(LonO + offLat))
    # And why does he add a null string to the end?

    f.write(datestr+' C K'+IonRadiusstr+' '+Latstr+' '+Lonstr+'')
    f.close()
    origdir = os.getcwd()
    os.chdir(IGRFpath)
    #os.system(IGRFpath+'geomag70.exe '+IGRFpath+'IGRF11.COF f '+IGRFpath+'input.txt '+str(path)+'ionFR/IGRF/geomag70_linux/output.txt')
    os.system('./geomag70.exe IGRF11.COF f input.txt output.txt')
    #g = open(''+str(path)+'ionFR/IGRF/geomag70_linux/output.txt', 'r')
    g = open('output.txt', 'r')
    data = g.readlines()
    print data
    g.close()
    #os.system('rm '+str(path)+'ionFR/IGRF/geomag70_linux/input.txt; rm '+str(path)+'ionFR/IGRF/geomag70_linux/output.txt')
    os.chdir(origdir)
    Xfield = abs(float(data[1].split()[10]))
    Yfield = abs(float(data[1].split()[11]))
    Zfield = abs(float(data[1].split()[12]))
    Xfield = Xfield*pow(10,-9)*Tesla2Gauss
    Yfield = Yfield*pow(10,-9)*Tesla2Gauss
    Zfield = Zfield*pow(10,-9)*Tesla2Gauss
    Totfield = Zfield*np.cos(ZenPunct) + Yfield*np.sin(ZenPunct)*np.sin(AzPunct) + Xfield*np.sin(ZenPunct)*np.cos(AzPunct)

    return Totfield
		






		
		





