from __future__ import print_function
import os
import sys
import subprocess
import numpy as np
import healpy as hp
from astropy import units as u
from astropy import constants as c
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle, Latitude, Longitude
import radiono as rad

if __name__ == '__main__':
    # PAPER INFO
    nside = 16
    npix = hp.nside2npix(nside)
    ipix = np.arange(npix)
    theta, phi = hp.pix2ang(nside, ipix)

    ra_strs = ('16h50m04.0s',)
    dec_strs = ('+79d11m25.0s',)
    lat_str = '52d54m54.64sn'
    lon_str = '6d36m16.04se'
    time_strs = ('2004-05-19T00:00:00',)
    height = 0

    #

    for time_str in time_strs:
        for ra_str, dec_str in zip(ra_strs, dec_strs):
            RM_dir = os.path.join(rad.base_path, 'RM_files/{date}'.format(date=time_str.split('T')[0]))
            if not os.path.exists(RM_dir):
                os.makedirs(RM_dir)

            year, month, day = time_str.split('T')[0].split('-')
            IONEX_file = rad.IONEX_file_needed(year, month, day)
            IONEX_name = os.path.join(rad.base_path, IONEX_file)

            lat_obs = Latitude(Angle(lat_str[:-1]))
            lon_obs = Longitude(Angle(lon_str[:-1]))

            start_time = Time(time_str)

            location = EarthLocation(lat=lat_obs, lon=lon_obs, height=height * u.m)

            TEC, _, all_info = rad.read_IONEX_TEC(IONEX_name)

            a, rms_a, ion_height = all_info[7:]

            tec_hp = rad.interp_time(a, TEC['lat'], TEC['lon'])
            rms_hp = rad.interp_time(rms_a, TEC['lat'], TEC['lon'])

            # predict the ionospheric RM for every hour within a day 
            UTs = np.linspace(0, 23, num=24)

            for UT in UTs:
                hour = rad.std_hour(UT)
                ra_dec = SkyCoord(ra=ra_str, dec=dec_str,
                                  location=location, obstime=start_time + UT * u.hr)
                altaz = ra_dec.altaz

                alt_src = altaz.alt
                az_src = altaz.az
                zen_src = altaz.zen

                off_lat, off_lon,\
                az_punct, zen_punct = rad.punct_ion_offset(lat_obs.radian,
                                                           az_src.radian,
                                                           zen_src.to(u.radian).value,
                                                           ion_height)
                coord_lat, coord_lon = rad.get_coords(lat_str, lon_str,
                                                      lat_obs, lon_obs,
                                                      np.degrees(off_lat),
                                                      np.degrees(off_lon))
                B_para = rad.B_IGRF(year, month, day,
                                    coord_lat, coord_lon,
                                    ion_height, az_punct, zen_punct)

                TEC_path, RMS_TEC_path = rad.interp_space(tec_hp[UT], rms_hp[UT],
                                                          coord_lat, coord_lon,
                                                          zen_punct)

                new_file = os.path.join(RM_dir, 'IonRM{hour}.txt'.format(hour=hour))
                rad.get_results(hour, new_file, B_para, TEC_path, RMS_TEC_path)
