'''
radiono

purpose | Module used to gather information from IONEX files

Functions
---------
std_hour | converts hour into consistent string representation
ion_RM | writes ionospheric RM to files for future use
'''
from __future__ import print_function
import os
import healpy as hp
import numpy as np
from astropy import constants as c
import astropy.coordinates as coord
from astropy import units
from astropy.time import Time
from radiono import physics as phys, interp as itp, ionex_file as inx

rad_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.abspath(os.path.join(rad_path, '..'))
rm_dir = os.path.join(root_dir, 'RM_files')
ionex_dir = os.path.join(root_dir, 'TEC')
TECU = 1e16
TEC2m2 = 0.1 * TECU
earth_radius = c.R_earth.value #6371000.0 # in meters
tesla_to_gauss = 1e4

def get_rm_map(date_str, verbose=False):
    """
    Returns the healpix RM map at resolution nside=16 in coordinates of RA/dec
    for the specified date at the HERA/PAPER site:
    (lat='-30d43m17.5s', lon='21d25m41.9s')

    Parameters:
        date_str: string denoting the date, order year-month-day e.g. date_str = '2012-06-01'

    Returns:
        RM_map: numpy array with RM_map.shape = (24, 3072), where each RM_map[i,:]
            is a healpix map of RMs as at the HERA/PAPER site, in J2000 RA/DEC coordinates.
    """

    def RM(B_para, TEC_path):
        IFR = 2.6e-17 * B_para * TEC_path
        return IFR

    filedir = os.path.join(base_path, 'RM_maps')
    filepath = os.path.join(filedir, date_str + '.npz')

    if os.path.exists(filepath):
        if verbose:
            print('Restoring ' + date_str + ' from save file.')
        RM_map = np.load(filepath)['RM_map']
        return RM_map

    year, month, day = date_str.split('-')

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

    lat_str = str(c_icrs.location.latitude)[1:] + 's'
    lon_str = str(c_icrs.location.longitude) + 'e'

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

        TEC_path[t], _ = itp.interp_space(tec_hp[t], rms_hp[t], lat, lon, za_p)

    RM_map = RM(B_para, TEC_path)

    if verbose:
        print('Saving ' + date_str + '.')
    np.savez(filepath, RM_map=RM_map)

    return RM_map

if __name__ == '__main__':
    print('This is not a script anymore')
