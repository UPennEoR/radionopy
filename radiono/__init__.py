'''
radiono

purpose | Module used to gather information from IONEX files

Functions
---------
std_hour | converts hour into consistent string representation
write_RM | writes ionospheric RM to files for future use
write_radec | writes RAs and DECs to file
maps_to_npz | writes maps to npz files
'''
from __future__ import print_function
import os
import numpy as np
from astropy import constants as c, units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

rad_dir = os.path.dirname(__file__)
root_dir = os.path.abspath(os.path.join(rad_dir, '..'))
rm_dir = os.path.join(root_dir, 'RM_files')
ionex_dir = os.path.join(root_dir, 'TEC')
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

def write_RM(hour, new_file, B_para, TEC_path, RMS_TEC_path, write_to_file=True):
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

def write_radec(UT, radec_file, alt_src, az_src, date_str, location, height=1051, verbose=False):
    '''
    writes ra, dec to file

    Parameters
    ----------
    UT | int: numbered hour of time
    radec_file | str: name of file to write ra, dec to
    alt_src | array: array of altitudes
    az_src | array: array of azimuths
    location | object: location object
    height | Optional[float]: height observation taken at
    verbose | Optional[bool]: whether to print values or not
    '''
    hour = std_hour(UT, verbose=False)

    start_time = Time(date_str)

    altaz = SkyCoord(alt=alt_src * u.deg, az=az_src * u.deg,
                     location=location, obstime=start_time + UT * u.hr,
                     frame='altaz')
    ra = altaz.icrs.ra
    dec = altaz.icrs.dec

    with open(radec_file, 'w') as f:
        for r, d in zip(ra, dec):
            f.write('{ra} {dec}\n'.format(ra=r.value, dec=d.value))

def maps_to_npz(time_str, npix, loc_str='PAPER', verbose=True):
    '''
    writes maps to npz files

    Parameters
    ----------
    time_str | str: time in string representation
    npix | int: number of pixels for healpix map
    loc_str | str: name of location to incorporate into filename output
    verbose | Optional[bool]: whether to print values or not
    '''
    #I could fish around in the file read to get npix and then re-loop, but why not just be lazy sometimes
    rng = np.arange(24)
    final_TEC, final_rm, final_drm, ra, dec = np.zeros((rng.shape[0], npix)),\
                                              np.zeros((rng.shape[0], npix)),\
                                              np.zeros((rng.shape[0], npix)),\
                                              np.zeros((rng.shape[0], npix)),\
                                              np.zeros((rng.shape[0], npix))
    RM_dir = os.path.join(rm_dir, '{date}'.format(date=time_str.split('T')[0]))
    for UT in rng:
        rm_file = os.path.join(RM_dir, 'IonRM{num}.txt'.format(num=std_hour(UT, verbose=verbose)))
        radec_file = os.path.join(RM_dir, 'radec{num}.txt'.format(num=std_hour(UT, verbose=verbose)))
        
        _, TEC, B, RM, dRM = np.loadtxt(rm_file, unpack=True)
        RA, DEC = np.loadtxt(radec_file, unpack=True)
        
        final_TEC[UT, :] = TEC
        final_rm[UT, :] = RM
        final_drm[UT, :] = dRM
        ra[UT,:] = RA
        dec[UT,:] = DEC

    f_name = ''.join((time_str.split('T')[0], '_', loc_str, '.npz'))
    npz_dir = os.path.join(root_dir, 'npz')
    if not os.path.exists(npz_dir):
        os.mkdir(npz_dir)
    npz_file = os.path.join(npz_dir, f_name)
    if verbose:
        print('Saving TEC, RM +/- dRM data and RA/Dec mapping to {filename}'.format(filename=npz_file))

    np.savez(npz_file, TEC=final_TEC, RM=final_rm, dRM=final_drm, RA=ra, DEC=dec)

if __name__ == '__main__':
    print('This is not a script anymore')
