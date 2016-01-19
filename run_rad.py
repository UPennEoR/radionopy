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

 # Defining some variables for further use
TECU = pow(10, 16)
TEC2m2 = 0.1 * TECU
earth_radius = c.R_earth.value #6371000.0 # in meters <-- if we use the hard-coded version, we replicate ionFR results to high precision
tesla_to_gauss = pow(10, 4)

o = optparse.OptionParser()
o.set_usage('do_radionopy.py [options] output_filename')
o.set_description(__doc__)
o.add_option('--basepath',dest='basepath',type='string',help='Base path for IONEX file and IONEX fortran executables. Default is current working directory. If you give path in the -f option, be sure to put a blank string in here.') #XXX NEED to improve this functionality

o.add_option('-f','--file',dest='IONEX_file',type='string',help='IONEX file we are getting TEC values from')
o.add_option('-t','--time',dest='time_str',type='string',help='report time in format "YYYY-MM-DDTHH:MM:SS" (adheres to IONEX format; needs to match the date the IONEX file correponds to)')
o.add_option('-n','--nside',dest='nside',type='int', help='healpix nside for maps',default=128)
o.add_option('--lon',dest='lon_str',type='string',help='longitude of observatory in form "DDdMMmSS.SSse (or sw)"',default='21d25m41.9se')
o.add_option('--lat',dest='lat_str',type='string',help='latitude of observatory in same form as longitude (sn or ss)',default='30d43m17.5ss')

o.add_option('-v','--verbose',default=False,action='store_true', help='Turn on more print statements')

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

_alt,_az,_altaz = [],[],[]
print 'Calculating Alt/Az per hour'
for a in altaz:
    _temp = a.transform_to(AltAz(obstime=start_time+UTs*u.hour,location=location))
    _alt.append(_temp.alt) #24 numbers*u.deg 
    _az.append(_temp.az)#24 numbers*u.deg
    _altaz.append(_temp.altaz) #1 SkyCoord object (with all metadata PER LOOP)

import IPython; IPython.embed()

#XXX inefficient to do same call twice. Is it easy to change read_IONEX_TEC to return both at once? If not, OK.
TEC, info = read_IONEX_TEC(IONEX_name)
RMS_TEC, rms_info = read_IONEX_TEC(IONEX_name, rms=True)

# Reading the altitude of the Ionosphere in km (from IONEX file)
hgt_ion = TEC['AltIon']


for i, UT in enumerate(UTs):
    if opts.v: print(UT)
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


"""
def get_ionRM_dict(lat_obs,lon_obs,alt_src,az_src,height_ion,time):
    #XXX how do I get the time from the altaz? 
    UTs = np.array(range(24))
    
    if alt_source.degree < 0: 



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
"""