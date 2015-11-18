import numpy as np
import pylab as plt
import healpy as hp
from astropy import units as u
from astropy import constants as c
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle, Latitude, Longitude

import radionopyBeta as ri
reload(ri)

# Nominally try to reproduce the output of this command
# ionFRM.py 16h50m04.0s+79d11m25.0s 52d54m54.64sn 6d36m16.04se 2004-05-19T00:00:00 CODG1400.04I
# Echo back what he has ... 
RAstr = '16h50m04.0s'
DECstr = '+79d11m25.0s'
LatStr = '52d54m54.64sn'
LonStr = '6d36m16.04se'
TimeStr = '2004-05-19T00:00:00' # This will actually work as input to the astropy Time function
IONEXfile = 'CODG1400.04I'

# Just do a dumb read and see if you can parse it ...
myTEC = ri.readIONEX(IONEXfile)

# ---------------------
# For future use
year,month,day = (TimeStr.split('T')[0]).split('-')

# For myself, I'm not quite sure what the right way to handle the
# recording of Lat/Lon - the annoying thing is parsing E/W, N/S for
# signs
LatObs = Latitude(Angle(LatStr[:-1]))
LonObs = Longitude(Angle(LonStr[:-1]))

# create the astropy location object to convert given RA/Dec to 
location = EarthLocation(lat=LatObs,lon=LonObs,height=0*u.m)
time0 = Time(TimeStr)

# Create a sky coordinate object, from which we can subsequently derive the necessary alt/az
radec = SkyCoord(ra=RAstr,dec=DECstr,location=location,obstime=time0)

# Read the TEC file once.  There is apparently only one ionosphere
# altitude per day (not per hour: CHECK.  Certainly the DHGT keyword
# only appears once)
TEC = ri.readIonexTEC(IONEXfile)
# Ionosphere altitude comes from the IONEX file
AltIon = TEC['AltIon']
# Now loop through times to reproduce ionFRM.py output
f = open(''+str(os.getcwd())+'/IonRM.txt', 'w')
UTs = np.linspace(0,23,num=24)
#for i,UT in enumerate(UTs):    
#    # Adjust the time to the current UT
#    radec = SkyCoord(ra=RAstr,dec=DECstr,location=location,obstime=time0+UT*u.hr)
#    # Calculate alt and az
#    AltSource = radec.altaz.alt
#    AzSource = radec.altaz.az
#    # ZenSource is a different kind o' object than Alt/Az
#    ZenSource = radec.altaz.zen
#    if (AltSource.degree > 0):
#        print i,AltSource,AzSource
#        # Calculate the ionospheric piercing point.  Inputs and outputs in radians
#        offLat,offLon,AzPunct,ZenPunct = ri.PuncIonOffset(LatObs.radian,AzSource.radian,ZenSource.to(u.radian).value,AltIon)
#        print offLat,offLon,AzPunct,ZenPunct
#        # Now we need to calculate the TEC at the given time at the
#        # requested Lat/Lon plus offsets.  This original code is
#        # goddamn mess.
#        VTEC = ri.interpTEC(TEC,LatO + offLat,LonO + offLon,UT)
#        TECpath = VTEC*TEC2m2/np.cos(ZenPunct)
#        VRMSTEC
#        RMSTECpath = VRMSTEC*TEC2m2/math.cos(ZenPunct)
#
#        # Get the mangetic field
#        Totfield = ri.B_IGRF(year,month,day,AltIon,LonObs,offLon,LatObs,offLat,AzPunct,ZenPunct)
#
#        # Saving the Ionosheric RM and its corresponding
#        # rms value to a file for the given 'hour' value
#        IFR = 2.6*pow(10,-17)*Totfield*TECpath
#        RMSIFR = 2.6*pow(10,-17)*Totfield*RMSTECpath
#        f.write(''+str(hour)+' '+str(TECpath)+' '+str(Totfield)+' '+str(IFR)+' '+str(RMSIFR)+'\n')
#
#f.close()




