'''
radiono.interp

authors | James Aguirre, Immanuel Washington, Saul Kohn

purpose | Module used to gather information from IONEX files

Functions
---------
interp_hp_time | interpolates healpix map in time
interp_time | converts square map into healpix map and interpolates in time
interp_space | interpolates healpix map in space
healpixellize | convert square map into healpix map
rotate_healpix_map | rotates healpix map
'''
from __future__ import print_function
import numpy as np
import healpy as hp
import radiono as rad

def interp_hp_time(map_i, map_j, t_i, t_j, t):
    '''
    interpolated healpix map in time

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
        print(t_i, t, t_j)
        print('Times will not work')
        return None

    w_i = float(t_j - t) / (t_j - t_i)
    w_j = float(t - t_i) / (t_j - t_i)
    dt_i_deg = -np.abs((t - t_i) * 360. / 24.)
    dt_j_deg = np.abs((t - t_j) * 360. / 24.)

    interp_map = w_i * rotate_healpix_map(map_i, [dt_i_deg, 0]) +\
                 w_j * rotate_healpix_map(map_j, [dt_j_deg, 0])

    return interp_map

def interp_time(maps, lat, lon, verbose=True):
    '''
    convert square map into healpix map
    interpolate healpix map in time

    Parameters
    ----------
    maps | array: square TEC maps
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
    nside = 16

    map_len = (len(maps) - 1) * 2
    hp_maps = []
    even_maps = [healpixellize(sq_map, lat_rad, lon_rad, nside, verbose=verbose)\
                 for sq_map in maps]

    for i, even_map in enumerate(even_maps[:-1]):
        hp_maps.append(even_map)
        if i < len(even_maps[:-1]):
            start_time = i
            end_time = (i + 2) % map_len
            mid_time = (end_time + start_time) / 2

            odd_map = interp_hp_time(even_maps[i], even_maps[i + 1],
                                     start_time, end_time, mid_time)
            hp_maps.append(odd_map)

    return np.array(hp_maps)

def interp_space(tec_hp, rms_hp, coord_lat, coord_lon, zen_punct):
    '''
    interpolates healpix maps in space

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

    TEC_path = np.array(VTEC) * rad.TEC2m2 / np.cos(zen_punct) # from vertical TEC to line of sight TEC
    RMS_TEC_path = np.array(VRMS_TEC) * rad.TEC2m2 / np.cos(zen_punct) # from vertical RMS_TEC to line of sight RMS_TEC

    return TEC_path, RMS_TEC_path

def healpixellize(f_in, theta_in, phi_in, nside, verbose=False):
    '''
    A dumb method for converting data f sampled at points theta and phi
    (not on a healpix grid)
    into a healpix at resolution nside

    Parameters
    ----------
    f_in | array: square map
    theta_in | array: latitude to grid to
    phi_in | array in: longitude to grid to
    nside | int: healpix value
    verbose | Optional[bool]: whether to print values or not

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

if __name__ == '__main__':
    print('This is not a script anymore')
