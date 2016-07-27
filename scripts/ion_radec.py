'''
radiono.scripts.ion_radec

authors | James Aguirre, Immanuel Washington, Saul Kohn

purpose | script to generate RM data from IONEX file using RAs and DECs
'''
from __future__ import print_function
import os
import numpy as np
import healpy as hp
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, Angle, Latitude, Longitude
import radiono as rad
from radiono import physics as phys, interp as itp, ionex_file as inx

def ion_RM(ra_strs, dec_strs, lat_str, lon_str, time_strs, height=0, ionex_dir=None, rm_dir=None):
    '''
    Parameters
    ----------
    ra_strs | list[str]: list of RAs for observation, corresponds by element to dec_strs
    dec_strs | list[str]: list of DECs for observation, corresponds by element for ra_strs
    lat_str | str: latitude that the observation was taken at
    lon_str | str: longitude that the observation was taken at
    time_strs | list[str]: list of dates for observation
    height | Optional[float]: height observation taken at in meters
    ionex_dir | str: directory in which ionex files are / will be located
    rm_dir | str: directory in which RM data files are / will be located
    '''
    lat_obs = Latitude(Angle(lat_str[:-1]))
    lon_obs = Longitude(Angle(lon_str[:-1]))
    location = EarthLocation(lat=lat_obs, lon=lon_obs, height=height * u.m)

    for time_str in time_strs:
        start_time = Time(time_str)
        RM_dir = os.path.join(rm_dir, '{date}'.format(date=time_str.split('T')[0]))
        if not os.path.exists(RM_dir):
            os.makedirs(RM_dir)
        for ra_str, dec_str in zip(ra_strs, dec_strs):
            year, month, day = time_str.split('T')[0].split('-')
            IONEX_file = inx.IONEX_file_needed(year, month, day)

            TEC, _, all_info = inx.read_IONEX_TEC(IONEX_file)

            a, rms_a, ion_height = all_info[7:]

            tec_hp = itp.interp_time(a, TEC['lat'], TEC['lon'])
            rms_hp = itp.interp_time(rms_a, TEC['lat'], TEC['lon'])

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

                if len(alt_src.shape) == 0:
                    alt_src = np.array([alt_src.item().degree])
                    az_src = np.array([az_src.item().degree])
                    zen_src = np.array([zen_src.value])

                coord_lat, coord_lon,\
                az_punct, zen_punct = phys.ipp(lat_str, lon_str,
                                               az_src, zen_src, ion_height)

                B_para = phys.B_IGRF(year, month, day,
                                     coord_lat, coord_lon,
                                     ion_height, az_punct, zen_punct)

                TEC_path, RMS_TEC_path = itp.interp_space(tec_hp[UT], rms_hp[UT],
                                                          coord_lat, coord_lon,
                                                          zen_punct)

                new_file = os.path.join(RM_dir, 'IonRM{hour}.txt'.format(hour=hour))
                rad.ion_RM(hour, new_file, B_para, TEC_path, RMS_TEC_path)

if __name__ == '__main__':
    ra_strs = ('16h50m04.0s',)
    dec_strs = ('+79d11m25.0s',)
    lat_str = '52d54m54.64sn'
    lon_str = '6d36m16.04se'
    time_strs = ('2004-05-19T00:00:00',)
    height = 0
    rm_dir = os.path.join(rad.base_path, 'RM_files')

    ion_RM(ra_strs, dec_strs, lat_str, lon_str, time_strs, height=height, ionex_dir=rad.ionex_dir, rm_dir=rm_dir)
