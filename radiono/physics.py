'''
radiono.physics

purpose | Module used to gather information from IONEX files

Functions
---------
B_IGRF | uses geomag.c script to acquire the B field for specific day
punct_ion_offset | finds ionosphere offsets and puncture points
get_coords | converts coordinates to degrees and includes offsets
ipp | acquires correct coordinates and puncture points
'''
from __future__ import print_function
import os
import datetime
import subprocess
import numpy as np
import healpy as hp
from astropy.coordinates import Angle, Latitude, Longitude
from astropy import constants as c
import radiono as rad
import tempfile

from random import choice
from string import ascii_lowercase

earth_radius = c.R_earth.value #6371000.0 # in meters
tesla_to_gauss = 1e4

def B_IGRF(year, month, day, coord_lat, coord_lon, ion_height, az_punct, zen_punct, mag_dir='IGRF/geomag70_linux', mag_file='IGRF11.COF'):
    '''
    calculates the B field for a particular date at particular coordinates

    Parameters
    ----------
    year | int: year
    month | int: numbered month of the year
    day | int: numbered day of the month
    coord_lat | array[float]: array of latitudes
    coord_lon | array[float]: array of longitudes
    ion_height | float: ionosphere height
    az_punct | array: array of azimuth puncture points
    zen_punct | array: array of zenith puncture points

    Returns
    -------
    array: B field at each coordinate
    '''

    script_dir = os.path.join(rad.root_dir, 'IGRF/geomag70_linux/')

    def random_io_file(script_dir):
        """
        Create a file with a unique name for use by the IGRF script.
        """

        # generate a random string of lowercase letters a-z
        io_name = ''.join(choice(ascii_lowercase) for i in range(10)) + '.txt'
        io_file = os.path.join(script_dir,io_name)

        # check to see if this file already exists, and generate names
        # until a new name is found
        while os.path.exists(io_file):
            io_name = ''.join(choice(ascii_lowercase) for i in range(10)) + '.txt'
            io_file = os.path.join(script_dir,io_name)

        # create an empty file with the new name
        open(io_file, 'w').close()

        return io_file, io_name

    input_file, input_name = random_io_file(script_dir)
    output_file, output_name = random_io_file(script_dir)

    # Calculation of the total magnetic field along the line of sight at the IPP
    sky_rad = (earth_radius + ion_height) / 1000.0
    with open(input_file, 'w') as f:
        for co_lat, co_lon in zip(coord_lat, coord_lon):
            f.write(('{year},{month},{day} '
                     'C K{sky_rad} '
                     '{ipp_lat} {ipp_lon}\n').format(year=year,
                                                     month=month,
                                                     day=day,
                                                     sky_rad=sky_rad,
                                                     ipp_lat=co_lat,
                                                     ipp_lon=co_lon))

    #runs the geomag exe script
    script_dir = os.path.join(rad.root_dir, 'IGRF/geomag70_linux/')
    working_dir = os.getcwd()
    os.chdir(script_dir)
    os.system('./geomag70 IGRF11.COF f ' + input_name + ' ' + output_name + ' > /dev/null 2>&1')
    os.chdir(working_dir)

    #record B_parallel field in numpy array
    B_para = []
    with open(output_file, 'r') as g:
        all_data = g.readlines()

        for i, data in enumerate(all_data[1:]):
            x_field,\
            y_field,\
            z_field = [float(field_data) * 1e-9 * tesla_to_gauss\
                       for field_data in data.split()[10:13]]
            #ABOVE CHANGE: abs(float(field_data)) -> float(field_data)
            B_paras = z_field * np.cos(zen_punct[i]) +\
                      y_field * np.sin(zen_punct[i]) * np.sin(az_punct[i]) +\
                      x_field * np.sin(zen_punct[i]) * np.cos(az_punct[i])

            B_para.append(B_paras)

    # clean-up the input and output files used by IGRF
    os.remove(input_file)
    os.remove(output_file)

    return np.array(B_para)

def punct_ion_offset(lat_obs, az_src, zen_src, ion_height):
    '''
    calculates offsets and puncture points

    Parameters
    ----------
    lat_obs | object: latitude
    az_src | array: array of azimuths in radians
    zen_src | array: array of zeniths in radians
    ion_height | float: ionosphere height

    Returns
    -------
    tuple:
        array: offset latitude
        array: offset longitude
        array: azimuth puncture array
        array: zenith puncture array
    '''
    #earth_radius = 6371000.0 # in meters

    # The 2-D sine rule gives the zenith angle at the
    # Ionospheric piercing point
    zen_punct = np.arcsin((earth_radius * np.sin(zen_src)) /\
                          (earth_radius + ion_height))

    # Use the sum of the internal angles of a triange to determine theta
    theta = zen_src - zen_punct

    # The cosine rule for spherical triangles gives us the latitude
    # at the IPP
    lat_ion = np.arcsin(np.sin(lat_obs) *\
                        np.cos(theta) + np.cos(lat_obs) *\
                        np.sin(theta) * np.cos(az_src))
    off_lat = lat_ion - lat_obs # latitude difference

    # Longitude difference using the 3-D sine rule (or for spherical triangles)
    off_lon = np.arcsin(np.sin(az_src) * np.sin(theta) / np.cos(lat_ion))

    # Azimuth at the IPP using the 3-D sine rule
    s_az_ion = np.sin(az_src) * np.cos(lat_obs) / np.cos(lat_ion)
    az_punct = np.arcsin(s_az_ion)

    return off_lat, off_lon, az_punct, zen_punct

def get_coords(lat_str, lon_str, lat_obs, lon_obs, off_lat, off_lon):
    '''
    converts coordinates into degrees
    adjusts for the offset

    Parameters
    ----------
    lat_str | str: latitude
    lon_str | str: longitude
    lat_obs | array: latitude objects
    lon_obs | array: longitude objects
    off_lat | array: offset latitude array
    off_lon | array: offset longitude array

    Returns
    -------
    tuple:
        array: array of offset fixed latitudes
        array: array of offset fixed longitudes
    '''
    if lat_str[-1] == 's':
        lat_val = -1
    elif lat_str[-1] == 'n':
        lat_val = 1
    if lon_str[-1] == 'e':
        lon_val = 1
    elif lon_str[-1] == 'w':
        lon_val = -1

    try:
        coord_lat = lat_val * (lat_obs.value + off_lat)
    except ValueError:
        print("The 'lat_str' variable is probably missing an 's' or 'n' suffix.")
    try:
        coord_lon = lon_val * (lon_obs.value + off_lon)
    except ValueError:
        print("Check strings.")

    return coord_lat, coord_lon

def ipp(lat_str, lon_str, az_src, zen_src, ion_height):
    '''
    acquires correct coordinates and puncture points

    Parameters
    ----------
    lat_str | str: latitude
    lon_str | str: longitude
    az_src | array: array of azimuths in degrees
    zen_src | array: array of zeniths in degrees
    ion_height | float: ionosphere height

    Returns
    -------
    tuple:
        array: array of offset fixed latitudes
        array: array of offset fixed longitudes
        array: azimuth puncture array
        array: zenith puncture array
    '''
    lat_obs = Latitude(Angle(lat_str[:-1]))
    lon_obs = Longitude(Angle(lon_str[:-1]))

    off_lat, off_lon, az_punct, zen_punct = punct_ion_offset(lat_obs.radian,
                                                             np.radians(az_src),
                                                             np.radians(zen_src),
                                                             ion_height)
    coord_lat, coord_lon = get_coords(lat_str, lon_str,
                                      lat_obs, lon_obs,
                                      np.degrees(off_lat), np.degrees(off_lon))

    return coord_lat, coord_lon, az_punct, zen_punct
