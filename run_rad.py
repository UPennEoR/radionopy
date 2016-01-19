"""
Script to execute radionopy functions, creating...
"""

import rad
import sys, os
import healpy as hp, numpy as np, optparse as o
from astropy import units as u
from astropy import constants as c
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle, Latitude, Longitude

 # Defining some variables for further use
TECU = pow(10, 16)
TEC2m2 = 0.1 * TECU
earth_radius = c.R_earth.value #6371000.0 # in meters <-- if we use the hard-coded version, we replicate ionFR results to high precision
tesla_to_gauss = pow(10, 4)

o = optparse.OptionParser()
o.set_usage('do_radionopy.py [options] output_filename')
o.set_description(__doc__)
o.add_option('--basepath',dest='base_path',type='string',help='Base path for IONEX file and IONEX fortran executables. Default is current working directory. If you give path in the -f option, be sure to put a blank string in here.') #XXX NEED to improve this functionality

o.add_option('-f','--file',dest='IONEX_file',type='string',help='IONEX file we are getting TEC values from')
o.add_option('-t','--time',dest='time_str',type='string',help='report time in format "YYYY-MM-DDTHH:MM:SS" (adheres to IONEX format; needs to match the date the IONEX file correponds to)')
o.add_option('-n','--nside',dest='nside',type='int', help='healpix nside for maps',default=128)

#XXX should change defaults to PAPER/HERA values (for warm fuzzy feeling)
o.add_option('--lon',dest='lon_str',type='string',help='longitude of observatory in form "DDdMMmSS.SSse (or sw)"',default='6d36m16.04se')
o.add_option('--lat',dest='lat_str',type='string',help='latitude of observatory in same form as longitude (sn or ss)',default='52d54m54.64sn')

o.add_option('-v','--verbose',help='Turn on more print statements')

opts,args = o.parse_args(sys.argv[1:])

#set output file
if len(args>1):
    print 'only using one (first) output_filename'
    args = args[0]

if opts.basepath != None: base_path = opts.basepath
else: base_path = os.getcwd()

with open(os.path.join(base_path, args), 'w') as f: pass #XXX why do we need the "with/as" statement?

print 'RADIONOPY'
print '    IONEX file:',opts.IONEX_file
print '    Time:',opts.time_str.split('T')
print '    nside:',opts.nside
print '    latitude:',opts.lat
print '    longitude:',opts.lon


IONEX_name = os.path.join(base_path, IONEX_file)
year, month, day = time_str.split('T')[0].split('-')
lon_obs = Longitude(Angle(lon_str[:-1]))
lat_obs = Latitude(Angle(lat_str[:-1]))
location = EarthLocation(lon=lon_obs, lat=lat_obs, height=0 * u.m) #XXX Karoo height?
start_time = Time(time_str)


### Here we need to accept an array of RA/Dec which correspond to the
### centers of healpix pixels, and all subsequent operations should
### allow ra/dec to be vectorized

### The returned value of the function is then an RA/Dec map of the RM
### above the array

### npix = hp.nside2npix(nside)
### ipix = np.arange(npix)
### ra,dec = hp.pix2ang(nside,ipix)
### ra_dec = SkyCoord() # go from healpix theta,phi radians to astropy ra,dec
### alt_source = ra_dec.altaz.al
### az_source = ra_dec.altaz.az
### that passes in to the new function


## Nominally try to reproduce the output of this command
## ionFRM.py 16h50m04.0s+79d11m25.0s 52d54m54.64sn 6d36m16.04se 2004-05-19T00:00:00 CODG1400.04I 
## Echo back what he has ... 

IONEX_name = os.path.join(base_path, IONEX_file)

year, month, day = time_str.split('T')[0].split('-')

lon_obs = Longitude(Angle(lon_str[:-1]))
lat_obs = Latitude(Angle(lat_str[:-1]))

location = EarthLocation(lon=lon_obs, lat=lat_obs, height=0 * u.m)
start_time = Time(time_str)

# Create a sky coordinate object, from which we can subsequently derive the necessary alt/az
ra_dec = SkyCoord(ra=ra_str, dec=dec_str, location=location, obstime=start_time)

TEC, info = read_IONEX_TEC(IONEX_name)
RMS_TEC, rms_info = read_IONEX_TEC(IONEX_name, rms=True)

# Reading the altitude of the Ionosphere in km (from IONEX file)
alt_ion = TEC['AltIon']

# predict the ionospheric RM for every hour within a day 
UTs = np.linspace(0, 23, num=24)
for i, UT in enumerate(UTs):
    print(UT)
    if UT < 10:
        hour = '0{hour}'.format(hour=int(UT))
    else:
        hour = '{hour}'.format(hour=int(UT))
    
    ra_dec = SkyCoord(ra=ra_str, dec=dec_str, location=location, obstime=start_time + UT * u.hr)

    # Calculate alt and az
    alt_source = ra_dec.altaz.alt
    az_source = ra_dec.altaz.az

    # zen_source is a different kind of object than Alt/Az
    zen_source = ra_dec.altaz.zen

### From here to the end should be a function that takes
### lat_obs,lon_obs,alt_src,az_src,height_ion and computes the RM,
### returning a dictionary with TEC, RMS_TEC, IFR, RMS_IFR,tot_field,
### etc ...
    if (alt_source.degree > 0):
        print(i, alt_source, az_source)
        # Calculate the ionospheric piercing point.  Inputs and outputs in radians
        off_lon, off_lat, az_punct, zen_punct = punct_ion_offset(lat_obs.radian, az_source.radian, zen_source.to(u.radian).value, alt_ion)
        print(off_lon, off_lat, az_punct, zen_punct)

        #coord_lon, coord_lat = get_coords(lon_str, lat_str, lon_obs, lat_obs, off_lon, off_lat)
        coord_lon, coord_lat = get_coords(lon_str, lat_str, lon_obs, lat_obs, off_lon * 180 / np.pi, off_lat * 180 / np.pi)

        TEC_path = TEC_paths(TEC, UT, coord_lon, coord_lat, zen_punct, info)
        RMS_TEC_path = TEC_paths(RMS_TEC, UT, coord_lon, coord_lat, zen_punct, rms_info)
        tot_field = B_IGRF(year, month, day, coord_lon, coord_lat, alt_ion, az_punct, zen_punct)

        # Saving the Ionosheric RM and its corresponding
        # rms value to a file for the given 'hour' value
        IFR = 2.6 * pow(10, -17) * tot_field * TEC_path
        RMS_IFR = 2.6 * pow(10, -17) * tot_field * RMS_TEC_path

        with open(os.path.join(base_path, 'IonRM.txt'), 'a') as f:
            f.write('{hour} {TEC_path} {tot_field} {IFR} {RMS_IFR}\n'.format(hour=hour, TEC_path=TEC_path, tot_field=tot_field,
                                                                             IFR=IFR, RMS_IFR=RMS_IFR))