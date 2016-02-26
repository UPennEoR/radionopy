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

# Defining some variables for further use
### Make the base path settable
base_path = os.path.expanduser('~/radionopy')
TECU = 1e16
TEC2m2 = 0.1 * TECU
earth_radius = c.R_earth.value #6371000.0 # in meters
tesla_to_gauss = 1e4

def B_IGRF(year, month, day, coord_lat, coord_lon, ion_height, az_punct, zen_punct):
    # Calculation of TEC path value for the indicated 'hour' and therefore 
    # at the IPP

    input_file = os.path.join(base_path, 'IGRF/geomag70_linux/input.txt')
    output_file = os.path.join(base_path, 'IGRF/geomag70_linux/output.txt')

    #uses lat_val, lon_val from above
    # Calculation of the total magnetic field along the line of sight at the IPP
    sky_rad = (earth_radius + ion_height) / 1000.0
    with open(input_file, 'w') as f:
        f.write(('{year},{month},{day} '
                 'C K{sky_rad} '
                 '{ipp_lat} {ipp_lon}\n').format(year=year,
                                                 month=month,
                                                 day=day,
                                                 sky_rad=sky_rad,
                                                 ipp_lat=coord_lat,
                                                 ipp_lon=coord_lon))

    #XXX runs the geomag exe script
    script_name = os.path.join('./', base_path, 'IGRF/geomag70_linux/geomag70')
    script_data = os.path.join(base_path, 'IGRF/geomag70_linux/IGRF11.COF')
    script_option = 'f'
    subprocess.call([script_name, script_data, script_option, input_file, output_file])

    with open(output_file, 'r') as g:
        all_data = g.readlines()

        x_field, y_field, z_field = [abs(float(field_data)) * 1e-9 * tesla_to_gauss for field_data in all_data[1].split()[10:13]]
        B_para = z_field * np.cos(zen_punct) +\
                  y_field * np.sin(zen_punct) * np.sin(az_punct) +\
                  x_field * np.sin(zen_punct) * np.cos(az_punct)

    return np.array(B_para)

def get_results(hour, TEC_path, RMS_TEC_path, B_para):
    # Saving the Ionosheric RM and its corresponding
    # rms value to a file for the given 'hour' value
    IFR = 2.6e-17 * B_para * TEC_path
    RMS_IFR = 2.6e-17 * B_para * RMS_TEC_path

    new_file = os.path.join(base_path, 'RM_files',
                                       'decra.txt')
    if hour == '00':
        try:
            os.remove(new_file)
        except:
            pass
    with open(new_file, 'a') as f:
        f.write(('{hour} {TEC_path} '
                 '{B_para} {IFR} '
                 '{RMS_IFR}\n').format(hour=hour,
                                       TEC_path=TEC_path,
                                       B_para=B_para,
                                       IFR=IFR,
                                       RMS_IFR=RMS_IFR))


if __name__ == '__main__':
    # PAPER INFO
    nside = 16
    npix = hp.nside2npix(nside)
    ipix = np.arange(npix)
    theta, phi = hp.pix2ang(nside, ipix)

    #alt_src = 90. - np.degrees(np.array(theta))
    #az_src = np.degrees(np.array(phi))
    #zen_src = np.degrees(np.array(theta))

    ra_str = '16h50m04.0s'
    dec_str = '+79d11m25.0s'
    lat_str = '52d54m54.64sn'
    lon_str = '6d36m16.04se'
    time_str = '2004-05-19T00:00:00' # This will actually work as input to the astropy Time function
    #lat_str = '30d43m17.5ss'
    #lon_str = '21d25m41.9se'
    #time_str = '2012-02-13T00:00:00'
    #IONEX_file = 'CODG0440.12I'
    height = 0 #1000

    #

    year, month, day = time_str.split('T')[0].split('-')
    IONEX_file = rad.IONEX_file_needed(year, month, day)
    IONEX_name = os.path.join(base_path, IONEX_file)

    lat_obs = Latitude(Angle(lat_str[:-1]))
    lon_obs = Longitude(Angle(lon_str[:-1]))

    start_time = Time(time_str)

    location = EarthLocation(lat=lat_obs, lon=lon_obs, height=height * u.m)

    TEC, _, all_info = rad.read_IONEX_TEC(IONEX_name)

    a, rms_a, ion_height = all_info[7:]

    tec_hp = rad.interp_time(a, TEC['lat'], TEC['lon'])
    rms_hp = rad.interp_time(rms_a, TEC['lat'], TEC['lon'])

    #off_lat, off_lon, az_punct, zen_punct = punct_ion_offset(lat_obs.radian,
    #                                                         np.radians(az_src),
    #                                                         np.radians(zen_src),
    #                                                         ion_height)
    #coord_lat, coord_lon = get_coords(lat_str, lon_str,
    #                                  lat_obs, lon_obs,
    #                                  np.degrees(off_lat), np.degrees(off_lon))

    #B_para = B_IGRF(year, month, day,
    #                coord_lat, coord_lon,
    #                ion_height, az_punct, zen_punct)

    # predict the ionospheric RM for every hour within a day 
    UTs = np.linspace(0, 23, num=24)

    for UT in UTs:
        hour = rad.std_hour(UT)    
        ra_dec = SkyCoord(ra=ra_str, dec=dec_str, location=location, obstime=start_time + UT * u.hr)
        altaz = ra_dec.altaz

        alt_src = altaz.alt
        az_src = altaz.az
        zen_src = altaz.zen

        off_lat, off_lon, az_punct, zen_punct = rad.punct_ion_offset(lat_obs.radian, az_src.radian, zen_src.to(u.radian).value, ion_height)
        coord_lat, coord_lon = rad.get_coords(lat_str, lon_str, lat_obs, lon_obs, np.degrees(off_lat), np.degrees(off_lon))
        B_para = B_IGRF(year, month, day,
                        coord_lat, coord_lon,
                        ion_height, az_punct, zen_punct)

        TEC_path, RMS_TEC_path = rad.interp_space(tec_hp[UT], rms_hp[UT],
                                                  coord_lat, coord_lon,
                                                  zen_punct)

        #results = {'TEC': TEC, 'RMS_TEC': RMS_TEC,
        #           'IFR': IFR, 'RMS_IFR': RMS_IFR,
        #           'B_para': B_para}
        get_results(hour, TEC_path, RMS_TEC_path, B_para)
