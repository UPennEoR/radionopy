'''
radiono.utils

purpose | helper functions for reading and writing IONEX/Healpix data

Functions
---------
std_hour | converts hour into consistent string representation
write_RM | writes ionospheric RM to file(s)
write_radec | writes RAs and DECs to file
'''
import ephem, numpy as np, healpy as hp
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

def eph2ionDate(date):
    """
    YYYY/MM/DD to YYYY-MM-DD
    """
    return '-'.join(date.split('/'))

def ion2ephDate(date):
    """
    YYYY-MM-DD to YYYY/MM/DD
    """
    return '/'.join(date.split('-'))

def ephemPAPER(date=None):
    """
    returns a ephem Observer object of the PAPER location
    optionally supply a date in string form 'YYYY/MM/DD'
    """
    site = ephem.Observer()
    site.lat,site.long,site.elevation = '-30.721527777777776','21.428305555555557',1000.
    if date is not None: site.date = date
    return site

def nextTransit(date,ra,dec,lat=-30.721527777777776,lon=21.428305555555557,elev=1000.):
    """
    Construct observer object (default PAPER site) and ask when the
    next transit of a given RA/Dec is.

    Output time is Universal Time (UT).
    Paramters
    ---------
    date | str: in form YYYY/MM/DD
    ra, dec, lat and lon | float: in degrees
    elevation | float: meters
    """
    #define where and when we are observing
    site = ephem.Observer()
    site.lat,site.long,site.elevation = str(lat),str(lon),elev
    site.date = date
    tp = ephem.FixedBody() # this is the point we are asking about
    tp._ra = np.radians(ra)
    tp._dec = np.radians(dec)
    tp._epoch = date
    tp.compute()
    tp_transit = site.next_transit(tp)
    return str(tp_transit)

def parseTransitBasic(trans_str,SunCheck=False):
    """
    This method can be used to find the location in the
    radionopy output array for a transit of a given
    pointing.

    Parameters
    ----------
    trans_str | str: output from nextTransit(...); string in form 'YYYY/MM/DD HH:MM:SS.ss'
    SunCheck | bool: should we check if the Sun is up? If it is, the returned tuple contains 'True'
    """
    _date,_time_UTC = trans_str.split()
    _date = '-'.join(_date.split('/'))
    _time = map(int,_time_UTC.split(':'))
    # zeroth-order estimation is to round to nearest UT
    if float(_time[1]) > 30: up = True
    else: up = False
    _hour = _time[0]
    if up: _hour+=1
    if not SunCheck: return (_date,_hour)
    else:
        sun = ephem.Sun()
        paper = ephemPAPER(date=trans_str)
        sun._epoch = trans_str
        sun.compute(paper)
        if float(repr(sun.alt)) > 0.: chk = True
        else: chk = False
        return (_date,_hour,chk,float(repr(sun.alt)))

def IndexToDeclRa(index,nside,deg=False):
    """
    Convert a HEALPix pixel to RA,Dec coordinates
    """
    theta,phi=hp.pixelfunc.pix2ang(nside,index)
    if deg: return -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)
    else: return theta-np.pi/2., 2.*np.pi-phi

def nsideToRaDec(nside):
    """
    Return two arrays (RA, Dec) based on a HEALPix
    grid of length 'nside'.
    """
    ipix = range(hp.nside2npix(nside))
    ras,decs = [],[]
    for p in ipix:
        dec,ra = IndexToDeclRa(p,nside)
        ras.append(ra)
        decs.append(dec)
    return np.array(ras),np.array(decs)
