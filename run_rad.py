"""
Script to execute radionopy functions, creating...
"""

import rad
import sys, os,  optparse
import healpy as hp, numpy as np
from astropy import units as u
from astropy import constants as c
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle, Latitude, Longitude

 

def get_results_dict(lat_str,lon_str,radec_src,height_ion,TECinfo,TEC_RMSinfo,verb=False):
    """
    Return the results of an radionopy run in a dictionary
    
    Requires: 
       lat_obs and lon_obs to be Astropy Latitude and Longitude objects
       radec_src to be an Astropy SkyCoord object
       height_ion to be a float
       TECinfo and RMSinfo XXX SHOULD NOT LIVE HERE
             
    XXX currently for an individual IPP, want to expand to give back healpix arrays
    XXX currently for an individual time stamp, could expand to incorporate the UT loop
    """
    earth_radius = c.R_earth.value #6371000.0 # in meters <-- if we use the hard-coded version, we replicate ionFR results to high precision
    
    results = {}
    kwds = ['TEC','RMS_TEC','RM','RMS_RM','Btot']
    
    #XXX recalculate lon_obs and lat_obs from strings? Need string form for the rad.get_coords method...
    lon_obs = Longitude(Angle(opts.lon_str[:-1]))
    lat_obs = Latitude(Angle(opts.lat_str[:-1]))
    # Calculate alt, az and zenith
    alt_src = radec_src.altaz.alt
    az_src = radec_src.altaz.az
    zen_src = radec_src.altaz.zen
    
    UT = radec_src.obstime.to_datetime() #a weakness that we inherit from ionFR is that the interpolation happens on integer hours. But this datetime object is flexible enough for our uses. 
    
    if alt_src.degree < 0: #below horizon -- give it nans
        for kw in kwds: results[kw] = np.nan
        return results
    
    if verb: print radec.obstime.value, alt_src, az_src
    # Calculate the ionospheric piercing point.  Inputs and outputs in radians
    off_lon, off_lat, az_punct, zen_punct = rad.punct_ion_offset(lat_obs.radian, az_src.radian, zen_src.to(u.radian).value, height_ion)
    
    if verb: print off_lon,off_lat,az_punct,zen_punct
    
    coord_lon, coord_lat = rad.get_coords(lon_str, lat_str, lon_obs, lat_obs, off_lon * 180 / np.pi, off_lat * 180 / np.pi)
    
    #XXX again, would be nice for TEC_paths to return both at once
    TEC_path, RMS_TEC_path = rad.TEC_paths(TEC, RMS_TEC, UT.hour, coord_lon, coord_lat, zen_punct, TECinfo. TEC_RMSinfo)
    
    tot_field = rad.B_IGRF(UT.year, UT.month, UT.day, coord_lon, coord_lat, height_ion, az_punct, zen_punct)
    
    RM = 2.6 * pow(10, -17) * tot_field * TEC_path
    RMS_RM = 2.6 * pow(10, -17) * tot_field * RMS_TEC_path
    
    results['TEC'] = TEC_path
    results['RMS_TEC'] = RMS_TEC_path
    results['Btot'] = tot_field
    results['RM'] = RM
    results['RMS_RM'] = RMS_RM
    return results

o = optparse.OptionParser()
o.set_usage('do_radionopy.py [options] output_filename')
o.set_description(__doc__)
o.add_option('--basepath',dest='basepath',type='string',help='Base path for IONEX file and IONEX fortran executables. Default is current working directory. If you give path in the -f option, be sure to put a blank string in here.') #XXX NEED to improve this functionality

o.add_option('-f','--file',dest='IONEX_file',type='string',help='IONEX file we are getting TEC values from')
o.add_option('-t','--time',dest='time_str',type='string',help='report time in format "YYYY-MM-DDTHH:MM:SS" (adheres to IONEX format; needs to match the date the IONEX file correponds to)')
o.add_option('-n','--nside',dest='nside',type='int', help='healpix nside for maps',default=128)

o.add_option('--lon',dest='lon_str',type='string',help='longitude of observatory in form "DDdMMmSS.SSse (or sw)"',default='21d25m41.9se')
o.add_option('--lat',dest='lat_str',type='string',help='latitude of observatory in same form as longitude (sn or ss)',default='30d43m17.5ss')

#RA and DEC options should be retired once we have full healpixellization working
o.add_option('--ra',dest='ra_str',type='string',help='R.A. of ionospheric piercing point in form HHhMMmSS.Ss',default='05h19m49.7s')
o.add_option('--dec',dest='dec_str',type='string',help='Declination of ionospheric piercing point in form -DDdMMmSS.Ss',default='-45d46m44s')

o.add_option('-v','--verbose',default=False,action='store_true', dest='verb', help='Turn on more print statements')

opts,args = o.parse_args(sys.argv[1:])

#set output file, parse path and output file options
if len(args)>1:
    print 'only using one (first) output_filename'
    args = args[0]
elif len(args)==0: raise ImplementationError('Need to supply output filename')
else: args = args[0]

if opts.basepath == None: base_path = os.getcwd()
else: base_path = opts.basepath

f = open(base_path+'/'+args, 'w')

print 'RADIONOPY'
print '    IONEX file:',opts.IONEX_file
print '    Time:',opts.time_str.split('T')
print '    nside:',opts.nside
print '    latitude:',opts.lat_str
print '    longitude:',opts.lon_str

#parse string options
IONEX_name = os.path.join(base_path, opts.IONEX_file)
year, month, day = opts.time_str.split('T')[0].split('-')
lon_obs = Longitude(Angle(opts.lon_str[:-1]))
lat_obs = Latitude(Angle(opts.lat_str[:-1]))
location = EarthLocation(lon=lon_obs, lat=lat_obs, height=0 * u.m) #XXX Karoo height?
start_time = Time(opts.time_str)

#start healpix funtimes
npix = hp.nside2npix(opts.nside)
ipix = np.arange(npix)
theta,phi = hp.pix2ang(opts.nside,ipix)
alt = (90. - np.degrees(np.array(theta))) * u.degree
az = (np.degrees(np.array(phi))) * u.degree

# Create a sky coordinate object which holds alt and az in degrees (along with some metadata) following healpix indexing
altaz = SkyCoord(alt=alt, az=az, obstime=start_time, frame='altaz', location=location)
radec = altaz.icrs

#split into hour, min and sec arrays
ra_h, ra_m, ra_s = radec.ra.hms
dec_h, dec_m, dec_s = radec.dec.dms


# predict the ionospheric RM for every hour within a day 
#UTs = np.linspace(0, 23, num=24)

UTs = np.array(range(24))
"""
_alt,_az,_altaz = [],[],[]
print 'Calculating Alt/Az per hour'
for a in altaz:
    _temp = a.transform_to(AltAz(obstime=start_time+UTs*u.hour,location=location))
    _alt.append(_temp.alt) #24 numbers*u.deg 
    _az.append(_temp.az)#24 numbers*u.deg
    _altaz.append(_temp.altaz) #1 SkyCoord object (with all metadata PER LOOP)
"""

#XXX inefficient to do same call twice. Is it easy to change read_IONEX_TEC to return both at once? If not, OK.
TEC, info = rad.read_IONEX_TEC(IONEX_name)
RMS_TEC, rms_info = rad.read_IONEX_TEC(IONEX_name, rms=True)

# Reading the altitude of the Ionosphere in km (from IONEX file)
hgt_ion = TEC['ion_height']


for i, UT in enumerate(UTs):
    if opts.verb: print(UT)
    if UT < 10:
        hour = '0{hour}'.format(hour=int(UT))
    else:
        hour = '{hour}'.format(hour=int(UT))
    
    #TESTING -- will replace with altaz.icrs arrays
    radec_src = SkyCoord(ra=opts.ra_str, dec=opts.dec_str, location=location, obstime=start_time + UT * u.hr)
    
    R = get_results_dict(opts.lat_str,opts.lon_str,radec_src,hgt_ion,info,rms_info,verb=opts.verb)
    
    #import IPython; IPython.embed()
