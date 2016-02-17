from __future__ import print_function
import os
import sys
import datetime
import ftplib
import shutil
import subprocess
import numpy as np
import pylab as plt
import healpy as hp
from astropy import units as u
from astropy import constants as c
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle, Latitude, Longitude
import rad
reload(rad)

# Defining some variables for further use
### Make the base path settable
base_path = os.path.expanduser('~/PyModules/radionopy')
TECU = 1e16
TEC2m2 = 0.1 * TECU
earth_radius = c.R_earth.value #6371000.0 # in meters
tesla_to_gauss = 1e4

nside = 16
npix = hp.nside2npix(nside)
ipix = np.arange(npix)
theta, phi = hp.pix2ang(nside, ipix)

alt = (90. - np.degrees(np.array(theta))) * u.degree
az = (np.degrees(np.array(phi))) * u.degree

lat_str = '30d43m17.5ss'
lon_str = '21d25m41.9se'
time_str = '2012-02-13T22:00:00'
#IONEX_file = 'CODG0440.12I'
height = 1000

#

year, month, day = time_str.split('T')[0].split('-')
IONEX_file = rad.IONEX_file_needed(year, month, day)
IONEX_name = os.path.join(base_path, IONEX_file)

lat_obs = -1.*Latitude(Angle(lat_str[:-1]))
lon_obs = Longitude(Angle(lon_str[:-1]))

start_time = Time(time_str)
location = EarthLocation(lat=lat_obs, lon=lon_obs, height=height * u.m)

TEC, RMS_TEC, all_info = rad.read_IONEX_TEC(IONEX_name)

info = all_info[:7] + (all_info[7],)
rms_info = all_info[:7] + (all_info[8],)
ion_height = all_info[9]

_, _, points_lat, _, _, points_lon, number_of_maps, a = info
_, _, _, _, _, _, _, rms_a = rms_info

newa = rad.interp_time(points_lat, points_lon, number_of_maps, 25, a)
rmsa = rad.interp_time(points_lat, points_lon, number_of_maps, 25, rms_a)

#TEC_t = rad.interp_time(points_lat, points_lon, number_of_maps, 25, a)
#interp_time(points_lat, points_lon, number_of_maps, total_maps, a):

alt_src = np.radians(alt.value)
az_src = np.radians(az.value)
zen_src = np.pi/2.-alt_src

off_lat, off_lon, az_punct, zen_punct = rad.punct_ion_offset(lat_obs.radian, az_src, zen_src, ion_height)

coord_lon = np.degrees(off_lon)+lon_obs.value
coord_lat = np.degrees(off_lat)+lat_obs.value

UTs = np.linspace(0, 23, num=24)

map0 = rad.TEC2HP(TEC,0)
map2 = rad.TEC2HP(TEC,1)

def interp_hp_time(map_i,map_j,t_i,t_j,t):
    # Need to check that
    if (not (t_i <= t <= t_j)):
        print("Nope")
        return
    w_i = (t_j - t)/(t_j - t_i)
    w_j = (t - t_i)/(t_j - t_i)
    dt_i_deg = -np.abs((t - t_i)*360./24.)
    dt_j_deg = np.abs((t - t_j)*360./24.)
    interp_map = w_i * hpt.rotate_healpix_map(map0,[dt_i_deg,0]) + w_j * hpt.rotate_healpix_map(map2,[dt_j_deg,0])
    return interp_map

interp_map = interp_hp_time(map0,map2,0.,2.,0.5)

hp.mollview(map0,title='Map 0')
hp.mollview(map2,title='Map 2')
hp.mollview(hpt.rotate_healpix_map(map0,[-15,0]),title='Map 0 +15')
hp.mollview(hpt.rotate_healpix_map(map2,[15,0]),title='Map 2 -15')
hp.mollview(interp_map,title='Interpolated')
plt.show()

#for UT in UTs:
#    hour = std_hour(UT)    




