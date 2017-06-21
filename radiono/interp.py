'''
radiono.interp

purpose | Module used to gather information from IONEX files

Functions
---------
interp_hp_time | interpolates healpix maps between times
ionex2healpix | converts square map into healpix map and interpolates in time to get ionosphere maps in healpix form and 1-hour time intervals
get_los_tec | finds the LoS TEC on a given healpix map
healpixellize | convert square map into healpix map
rotate_healpix_map | rotates healpix map
'''
from __future__ import print_function
import numpy as np
import healpy as hp
import radiono as rad
import utils
from bisect import bisect_right

TECU = 1e16
TEC2m2 = 0.1 * TECU

def interp_hp_time(map_i, map_j, t_i, t_j, t):
    '''
    interpolates healpix map between times

    Parameters
    ----------
    map_i | array: healpix map
    map_j | array: healpix map with different time than map_i
    t_i | float: time at which map_i is at
    t_j | float: time at which map_j is at
    t | time to interpolate map to which is between t_i and t_j

    Returns
    -------
    array: interpolated healpix map
    '''
    # Need to check that
    if not (t_i <= t <= t_j):
        raise ValueError('Time %f cannot be interpolated between %f and %f'%(t,t_i,t_j))

    w_i = float(t_j - t) / (t_j - t_i)
    w_j = float(t - t_i) / (t_j - t_i)
    dt_i_deg = -np.abs((t - t_i) * 360. / 24.)
    dt_j_deg = np.abs((t - t_j) * 360. / 24.)

    interp_map = w_i * rotate_healpix_map(map_i, [dt_i_deg, 0]) +\
                 w_j * rotate_healpix_map(map_j, [dt_j_deg, 0])

    return interp_map

def ionex2healpix(maps, UTs, lat, lon, verbose=False):
    '''
    convert square map into healpix map
    interpolate healpix map in time

    Parameters
    ----------
    maps | array: square TEC maps
    UTs | list of decmial hours at which to compute TEC maps
    lat | array[float]: array of latitudes
    lon | array[float]: array of longitudes
    verbose | Optional[bool]: whether to print values or not

    Returns
    -------
    array: time interpolated healpix map
    '''
    nlat = len(lat)
    nlon = len(lon)
    lat_rad = np.outer(np.radians(90. - lat), np.ones(nlon))
    lon_rad = np.outer(np.ones(nlat), np.radians(lon % 360))

    node_hours = range(0,26,2)
    node_intervals = list(zip(node_hours[:-1], node_hours[1:]))

    node_maps = [healpixellize(sq_map, lat_rad, lon_rad, verbose=verbose)\
                 for sq_map in maps]

    hp_maps = []
    for UT in UTs:
        i = bisect_right(node_hours, UT) - 1
        if (i in range(len(node_intervals))) is False:
            raise ValueError

        interval = node_intervals[i]

        s,f = interval

        map_UT = interp_hp_time(node_maps[i], node_maps[i+1], s, f, UT)
        hp_maps.append(map_UT)

    return np.array(hp_maps)


def get_los_tec(tec_hp, rms_hp, coord_lat, coord_lon, zen_punct):
    '''
    finds the LoS TEC on a given healpix map

    Parameters
    ----------
    tec_hp | array: tec healpix map
    rms_hp | array: rms tec healpix map
    coord_lat | array[float]: array of latitudes
    coord_lon | array[float]: array of longitudes
    zen_punct | array: array of zenith puncture points

    Returns
    -------
    tuple:
        array: space interpolated TEC map
        array: space interpolated RMS TEC map
    '''
    lat_rad = np.radians(90. - coord_lat)
    lon_rad = np.radians(coord_lon % 360)
    VTEC = hp.get_interp_val(tec_hp, lat_rad, lon_rad)
    VRMS_TEC = hp.get_interp_val(rms_hp, lat_rad, lon_rad)

    TEC_path = np.array(VTEC) * TEC2m2 / np.cos(zen_punct) # from vertical TEC to line of sight TEC
    RMS_TEC_path = np.array(VRMS_TEC) * TEC2m2 / np.cos(zen_punct) # from vertical RMS_TEC to line of sight RMS_TEC

    return TEC_path, RMS_TEC_path

def healpixellize(f_in, theta_in, phi_in, nside=16, verbose=False):
    '''
    A dumb method for converting data f sampled at points theta and phi
    (not on a healpix grid)
    into a healpix at resolution nside

    Parameters
    ----------
    f_in | array: square map
    theta_in | array: latitude to grid to
    phi_in | array in: longitude to grid to
    nside | Optional[int]: healpix value (IONEX maps give nside=16 natively)
    verbose | Optional[bool]: whether to print values and info or not

    Returns
    -------
    array: healpixellized version of square map input
    '''
    # Input arrays are likely to be rectangular, but this is inconvenient
    f = f_in.flatten()
    theta = theta_in.flatten()
    phi = phi_in.flatten()

    pix = hp.ang2pix(nside, theta, phi)

    hp_map = np.zeros(hp.nside2npix(nside))
    hits = np.zeros(hp.nside2npix(nside))

    for i, v in enumerate(f):
        # Find the nearest pixels to the pixel in question
        neighbours, weights = hp.get_interp_weights(nside, theta[i], phi[i])
        # Add weighted values to hp_map
        hp_map[neighbours] += v * weights
        # Keep track of weights
        hits[neighbours] += weights

    hp_map = hp_map / hits
    wh_no_hits = np.where(hits == 0)
    if verbose:
        print('pixels with no hits', wh_no_hits[0].shape)
    hp_map[wh_no_hits[0]] = hp.UNSEEN

    wh = np.where(hp_map == np.nan)[0]
    for i, w in enumerate(wh):
        neighbors = hp.get_interp_weights(nside, theta[i], phi[i])
        hp_map[w] = np.median(neighbors)

    return hp_map

def rotate_healpix_map(map_in, rot):
    '''
    Will rotate the pixels of a map into (effectively) a new ordering
    representing a rotation of the function.
    Not sure why this isn't implemented in healpy directly (maybe it is).
    In order to map each pixel exactly to a new one,
    the transform is only accurate to the pixel size.

    Parameters
    ----------
    map_in | array: healpix map
    rot | tuple[float, float]: rotation degrees

    Returns
    -------
    array: rotated healpix map
    '''
    npix = len(map_in)
    nside = hp.npix2nside(npix)

    rot_map = np.zeros(npix)
    ipix = np.arange(npix)
    theta, phi = hp.pix2ang(nside, ipix)

    rotator = hp.Rotator(rot=rot)

    # For each pixel in the new map, find where it would have come
    # from in the old
    theta_rot, phi_rot = rotator(theta, phi)
    ipix_rot = hp.ang2pix(nside, theta_rot, phi_rot)

    rot_map = map_in[ipix_rot]

    return rot_map
