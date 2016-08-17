'''
radiono.utils

purpose | helper functions for reading and writing IONEX/Healpix data

Functions
---------
std_hour | converts hour into consistent string representation
write_RM | writes ionospheric RM to file(s)
write_radec | writes RAs and DECs to file
'''

from astropy import constants as c, units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

def std_hour(UT, verbose=False):
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
