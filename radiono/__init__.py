'''
radiono

authors | James Aguirre, Immanuel Washington, Saul Kohn

purpose | Module used to gather information from IONEX files

Functions
---------
std_hour | converts hour into consistent string representation
ion_RM | writes ionospheric RM to files for future use
'''
from __future__ import print_function
import os
from astropy import constants as c

import astropy.coordinates as coord
from astropy import units
from astropy.time import Time

from . import physics as phys, interp as itp, ionex_file as inx

rad_path = os.path.dirname(os.path.realpath(__file__))
base_path = os.path.abspath(os.path.join(rad_path, '..'))
rm_dir = os.path.join(base_path, 'RM_files')
ionex_dir = os.path.join(base_path, 'TEC')
TECU = 1e16
TEC2m2 = 0.1 * TECU
earth_radius = c.R_earth.value #6371000.0 # in meters
tesla_to_gauss = 1e4

def std_hour(UT, verbose=True):
    '''
    converts hour into standardized string

    Parameters
    ----------
    UT | int: numbered hour of time
    verbose | Optional[bool]: whether to print values or not

    Returns
    -------
    str: standardized string hour
    '''
    if verbose:
        print(int(UT))
    if UT < 10:
        hour = '0{hour}'.format(hour=int(UT))
    else:
        hour = '{hour}'.format(hour=int(UT))

    return hour

def ion_RM(hour, new_file, B_para, TEC_path, RMS_TEC_path, write_to_file=True):
    '''
    writes ionospheric RM to file

    Parameters
    ----------
    hour | str: standardized hour
    new_file | str: filename to write RM to
    B_para | array: B field parallel to zenith
    TEC_path | array: line of sight TEC
    RMS_TEC_path | array: line of sight RMS TEC
    write_to_file | bool: whether or not to write to file
    '''
    # Saving the Ionosheric RM and its corresponding
    # rms value to a file for the given 'hour' value
    IFR = 2.6e-17 * B_para * TEC_path
    RMS_IFR = 2.6e-17 * B_para * RMS_TEC_path

    if write_to_file:
        with open(new_file, 'w') as f:
            for tp, tf, ifr, rms_ifr in zip(TEC_path, B_para, IFR, RMS_IFR):
                f.write(('{hour} {TEC_path} '
                         '{B_para} {IFR} '
                         '{RMS_IFR}\n').format(hour=hour,
                                               TEC_path=tp,
                                               B_para=tf,
                                               IFR=ifr,
                                               RMS_IFR=rms_ifr))
    return B_para, IFR, RMS_IFR

def get_rm_map(date_str):

    def RM(B_para, TEC_path):
        IFR = 2.6e-17 * B_para * TEC_path
        return IFR

    year, month, day = date_str.split('T')[0].split('-')

    tec_hp, rms_hp, ion_height = inx.IONEX_data(year, month, day, verbose=False)

    nside = 2**4
    npix = hp.nside2npix(nside)
    hpxidx = np.arange(npix)

    cza, ra = hp.pix2ang(nside, hpxidx)
    dec = np.pi/2. - cza

    c_icrs = coord.SkyCoord(ra=ra * units.radian, dec=dec * units.radian, frame='icrs')
    c_icrs.location = coord.EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s') # Location of HERA/PAPER
    c_icrs.obstime = Time(date_str + 'T12:00:00') # use the middle of the day as the local coordinates for all 24 hours. Close enough?
    c_altaz = c_icrs.transform_to('altaz')

    az = np.array(np.radians(c_altaz.az))
    alt = np.array(np.radians(c_altaz.alt))
    za = np.pi/2. - alt

    lat, lon, az_p, za_p = phys.ipp(lat_str, lon_str,
                                np.degrees(az), np.degrees(za),
                                ion_height)

    B_para = phys.B_IGRF(year, month, day,
                    lat, lon,
                    ion_height,
                    az_p, za_p)

    TEC_path = np.zeros((24,npix))
    # RMS_TEC_path = np.zeros((24,npix)) # just the TEC for now
    for t in range(0,24):
        hour = rad.std_hour(t, verbose=False)

        TEC_path[t], _ = itp.interp_space(tec_hp[t], rms_hp[t],
                                                 lat, lon,
                                                 za_p)

    return RM(B_para, TEC_path)

if __name__ == '__main__':
    print('This is not a script anymore')
