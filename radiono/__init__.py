'''
radionpy.radiono

authors | James Aguirre, Immanuel Washington, Saul Kohn

purpose | Module used to gather information from IONEX files

Functions
---------
IONEX_file_needed | finds correct IONEX file and uncompresses it if necessary
get_IONEX_file | downloads IONEX file from ftp server if not on local machine
gen_IOMEX_list | pulls specific info from parsed IONEX file
read_IONEX_TEC | parses  IONEX file
IONEX_data | gathers all relevant info from IONEX file for specific date
interp_hp_time | interpolates healpix map in time
interp_time | converts square map into healpix map and interpolates in time
interp_space | interpolates healpix map in space
healpixellize | convert square map into healpix map
B_IGRF | uses C script to acquire the B field for specific day
punct_ion_offset | finds ionosphere offsets and puncture points
get_coords | converts coordinates to degrees and includes offsets
ipp | acquires correct coordinates and puncture points
rotate_healpix_map | rotates healpix map
std_hour | converts hour into consistent string representation
get_results | writes ionospheric RM to files for future use
'''
from __future__ import print_function
import os
import sys
import datetime
import ftplib
import subprocess
import numpy as np
import healpy as hp
from astropy import units as u
from astropy import constants as c
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle, Latitude, Longitude

# Defining some variables for further use
### Make the base path settable
rad_path = os.path.dirname(os.path.realpath(__file__)) #XXX makes (fair) assumption that you're executing this in the radionopy directory
base_path = os.path.abspath(os.path.join(rad_path, '..'))
tec_dir = os.path.join(base_path, 'TEC')
TECU = 1e16
TEC2m2 = 0.1 * TECU
earth_radius = c.R_earth.value #6371000.0 # in meters
tesla_to_gauss = 1e4

#######################
### IONEX FUNCTIONS ###
#######################

def IONEX_file_needed(year, month, day):
    '''
    pulls correct IONEX file for date input
    decompresses if necessary

    Parameters
    ----------
    year | int: year
    month | int: numbered month of the year
    day | int: numbered day of the month

    Returns
    -------
    str: name of IONEX file requested
    '''
    time_str = '{year} {month} {day}'.format(year=year, month=month, day=day)
    day_of_year = datetime.datetime.strptime(time_str, '%Y %m %d').timetuple().tm_yday

    if day_of_year < 10:
        day_of_year = '00{day_of_year}'.format(day_of_year=day_of_year)
    elif 10 <= day_of_year < 100:
        day_of_year = '0{day_of_year}'.format(day_of_year=day_of_year)

    # Outputing the name of the IONEX file you require
    ionex_file = 'CODG{day_of_year}0.{year_end}I'.format(day_of_year=day_of_year, year_end=str(year)[2:4])
    ionex_file_z = ''.join((ionex_file, '.Z'))

    if not os.path.exists(os.path.join(tec_dir, ionex_file))\
    and not os.path.exists(os.path.join(tec_dir, ionex_file_z)):
        ionex_file_z = get_IONEX_file(ionex_file, year, month, day)
        subprocess.call(['uncompress', ionex_file_z])

    #sys.exit()
    return os.path.join(tec_dir, ionex_file)

def get_IONEX_file(IONEX_file, year, month, day):
    '''
    downloads IONEX file from ftp server

    Parameters
    ----------
    IONEX_file | str: name of IONEX file to be downloaded
    year | int: year
    month | int: numbered month of the year
    day | int: numbered day of the month

    Returns
    -------
    str: name of IONEX file pulled from ftp server

    '''
    server = 'ftp.unibe.ch'

    ftp_dir = os.path.join('aiub/CODE/', year)
    IONEX_file_Z = ''.join((IONEX_file, '.Z'))

    if not os.path.exists(tec_dir):
        os.mkdir(tec_dir)

    IONEX_file_X = os.path.join(tec_dir, IONEX_file_Z)

    getting_file_str = 'Retrieving {IONEX_file_X} for {day} {month} {year}'.format(IONEX_file_X=IONEX_file_X, day=day, month=month, year=year)
    print(getting_file_str)

    try:
        ftp = ftplib.FTP(server, 'anonymous', 'jaguirre@sas.upenn.edu')
        ftp.cwd(ftp_dir)
        ftp.retrbinary(' '.join(('RETR', IONEX_file_Z)), open(IONEX_file_X, 'wb').write)
        ftp.quit()
    except:
        print('No file exists?')
        os.remove(IONEX_file_X)

    return IONEX_file_X

def gen_IONEX_list(IONEX_list):
    '''
    pulls information from parsed IONEX file

    Parameters
    ----------
    IONEX_list | list[str]: parsed list of IONEX file data

    Returns
    -------
    tuple:
        list[str]: IONEX file data without header
        list: RMS IONEX file data without header
        int: number of maps
        float: ionospheric height
        float: first latitude
        float: last latitude
        float: latitude step
        float: first longitude
        float: last longitude
        float: step longitude
    '''
    add = False
    rms_add = False
    base_IONEX_list = []
    RMS_IONEX_list = []
    for file_data in IONEX_list[:-1]:
        if not file_data:
            continue
        if file_data.split()[-2:] == ['RMS', 'MAP']:
            add = False
            rms_add = True
        elif file_data.split()[-2:] == ['IN', 'FILE']:
            number_of_maps = float(file_data.split()[0])

        if file_data.split()[0] == 'END' and file_data.split()[2] == 'HEADER':
            add = True

        if rms_add:
            RMS_IONEX_list.append(file_data)
        if add:
            base_IONEX_list.append(file_data)

        if file_data.split()[-1] == 'DHGT':
            ion_h = float(file_data.split()[0])
        elif file_data.split()[-1] == 'DLAT':
            start_lat, end_lat, step_lat = [float(data_item) for data_item in file_data.split()[:3]]
        elif file_data.split()[-1] == 'DLON':
            start_lon, end_lon, step_lon = [float(data_item) for data_item in file_data.split()[:3]]

    return base_IONEX_list, RMS_IONEX_list, number_of_maps, ion_h, start_lat, end_lat, step_lat, start_lon, end_lon, step_lon

def read_IONEX_TEC(filename, verbose=True):
    '''

    Parameters
    ----------
    filename | str: name of IONEX file
    verbose | Optional[bool]: whether to print values or not

    Returns
    -------
    tuple:
        dict: TEC values
        dict: RMS TEC values
        tuple:
            float: first latitude
            float: last latitude
            int: amount of latitude points
            float: first longitude
            float: last longitude
            int: amount of longitude points
            int: number of maps
            array: TEC array
            array: RMS TEC array
            float: ionosphere height in meters
    '''
    # Reading and storing only the TEC values of 1 day
    # (13 maps) into a 3D array

    # Opening and reading the IONEX file into memory
    with open(filename, 'r') as read_file:
        linestring = read_file.read()
        IONEX_list = linestring.split('\n')

    # creating a new array without the header and only
    # with the TEC maps
    base_IONEX_list, RMS_IONEX_list, number_of_maps, ion_h,\
    start_lat, end_lat, step_lat,\
    start_lon, end_lon, step_lon = gen_IONEX_list(IONEX_list)

    # Variables that indicate the number of points in Lat. and Lon.
    points_lat = ((end_lat - start_lat) / step_lat) + 1
    points_lon = ((end_lon - start_lon) / step_lon) + 1
    if verbose:
        print(start_lat, end_lat, step_lat)
        print(start_lon, end_lon, step_lon)
        print(points_lat, points_lon)
    
    # What are the Lat/Lon coords?
    latitude = np.linspace(start_lat, end_lat, num=points_lat)
    longitude = np.linspace(start_lon, end_lon, num=points_lon)

    TEC_list = []
    # Selecting only the TEC values to store in the 3-D array
    for new_IONEX_list in (base_IONEX_list, RMS_IONEX_list):
        # 3D array that will contain TEC values only
        a = np.zeros((number_of_maps, points_lat, points_lon))

        counter_maps = 1
        for i in range(len(new_IONEX_list)):
            # Pointing to first map (out of 13 maps)
            # then by changing 'counter_maps' the other
            # maps are selected
            if new_IONEX_list[i].split()[0] == str(counter_maps) and new_IONEX_list[i].split()[-4] == 'START':
                # pointing the starting latitude
                # then by changing 'counter_lat' we select
                # TEC data at other latitudes within
                # the selected map
                counter_lat = 0
                new_start_lat = float(str(start_lat))
                for item_lat in range(int(points_lat)):
                    if new_IONEX_list[i + 2 + counter_lat].split()[0].split('-')[0] == str(new_start_lat)\
                    or '-' + new_IONEX_list[i + 2 + counter_lat].split()[0].split('-')[1] == str(new_start_lat):
                        # Adding to array 'a' a line of latitude TEC data
                        # we account for the TEC values at negative latitudes
                        counter_lon = 0
                        for count_num in range(3, 8):
                            list_index = i + count_num + counter_lat
                            for new_IONEX_item in new_IONEX_list[list_index].split():
                                a[counter_maps - 1, item_lat, counter_lon] = new_IONEX_item
                                counter_lon = counter_lon + 1
                    counter_lat = counter_lat + 6
                    new_start_lat = new_start_lat + step_lat
                counter_maps = counter_maps + 1

        TEC_list.append({'TEC': np.array(a), 'a': a})

    tec_a = TEC_list[0]['a']
    rms_a = TEC_list[1]['a']
    TEC =  {'TEC': TEC_list[0]['TEC'], 'lat': latitude, 'lon': longitude}
    RMS_TEC =  {'TEC': TEC_list[1]['TEC'], 'lat': latitude, 'lon': longitude}
    
    return TEC, RMS_TEC, (start_lat, step_lat, points_lat,\
                          start_lon, step_lon, points_lon,\
                          number_of_maps, tec_a, rms_a, ion_h * 1000.0)

def IONEX_data(year, month, day, verbose=True):
    '''
    gathers all relevant IONEX info from file for specific date

    Parameters
    ----------
    year | int: year
    month | int: numbered month of the year
    day | int: numbered day of the month
    verbose | Optional[bool]: whether to print values or not

    Returns
    -------
    tuple:
        array: tec healpix map
        array: rms tec healpix map
        float: ionosphere height in meters
    '''
    IONEX_file = IONEX_file_needed(year, month, day)
    IONEX_name = os.path.join(base_path, IONEX_file)
    TEC, _, all_info = read_IONEX_TEC(IONEX_name, verbose=verbose)

    tec_a, rms_a, ion_height = all_info[7:]

    tec_hp = interp_time(tec_a, TEC['lat'], TEC['lon'], verbose=verbose)
    rms_hp = interp_time(rms_a, TEC['lat'], TEC['lon'], verbose=verbose)

    return tec_hp, rms_hp, ion_height

###############################
### INTERPOLATION FUNCTIONS ###
###############################

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

    TEC_path = np.array(VTEC) * TEC2m2 / np.cos(zen_punct) # from vertical TEC to line of sight TEC
    RMS_TEC_path = np.array(VRMS_TEC) * TEC2m2 / np.cos(zen_punct) # from vertical RMS_TEC to line of sight RMS_TEC

    return TEC_path, RMS_TEC_path
    
def healpixellize(f_in, theta_in, phi_in, nside, verbose=True):
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

####################################
### PHYSICS/COORDINATE FUNCTIONS ###
####################################

def B_IGRF(year, month, day, coord_lat, coord_lon, ion_height, az_punct, zen_punct):
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
    # Calculation of TEC path value for the indicated 'hour' and therefore 
    # at the IPP

    input_file = os.path.join(base_path, 'IGRF/geomag70_linux/input.txt')
    output_file = os.path.join(base_path, 'IGRF/geomag70_linux/output.txt')

    #uses lat_val, lon_val from above
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

    #XXX runs the geomag exe script
    script_name = os.path.join('./', base_path, 'IGRF/geomag70_linux/geomag70')
    script_data = os.path.join(base_path, 'IGRF/geomag70_linux/IGRF11.COF')
    script_option = 'f'
    subprocess.call([script_name, script_data, script_option, input_file, output_file])

    B_para = []
    with open(output_file, 'r') as g:
        all_data = g.readlines()

        for i, data in enumerate(all_data[1:]):
            x_field,\
            y_field,\
            z_field = [abs(float(field_data)) * 1e-9 * tesla_to_gauss\
                       for field_data in data.split()[10:13]]

            B_paras = z_field * np.cos(zen_punct[i]) +\
                      y_field * np.sin(zen_punct[i]) * np.sin(az_punct[i]) +\
                      x_field * np.sin(zen_punct[i]) * np.cos(az_punct[i])

            B_para.append(B_paras)

    return np.array(B_para)

def punct_ion_offset(lat_obs, az_src, zen_src, ion_height):
    '''
    calculates offsets and puncture points

    Parameters
    ----------
    lat_obs | object: latitude
    az_src | array: array of azimuths
    zen_src | array: array of zeniths
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

    coord_lat = lat_val * (lat_obs.value + off_lat)
    coord_lon = lon_val * (lon_obs.value + off_lon)

    return coord_lat, coord_lon

def ipp(lat_str, lon_str, az_src, zen_src, ion_height):
    '''
    acquires correct coordinates and puncture points

    Parameters
    ----------
    lat_str | str: latitude
    lon_str | str: longitude
    az_src | array: array of azimuths
    zen_src | array: array of zeniths
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


######################################
### EXECUTION AND HELPER FUNCTIONS ###
######################################

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

def get_results(hour, new_file, B_para, TEC_path, RMS_TEC_path):
    '''
    writes ionospheric RM to file

    Parameters
    ----------
    hour | str: standardized hour
    new_file | str: filename to write RM to
    B_para | array: B field parallel to zenith
    TEC_path | array: line of sight TEC
    RMS_TEC_path | array: line of sight RMS TEC
    '''
    # Saving the Ionosheric RM and its corresponding
    # rms value to a file for the given 'hour' value
    IFR = 2.6e-17 * B_para * TEC_path
    RMS_IFR = 2.6e-17 * B_para * RMS_TEC_path

    with open(new_file, 'w') as f:
        for tp, tf, ifr, rms_ifr in zip(TEC_path, B_para, IFR, RMS_IFR):
            f.write(('{hour} {TEC_path} '
                     '{B_para} {IFR} '
                     '{RMS_IFR}\n').format(hour=hour,
                                           TEC_path=tp,
                                           B_para=tf,
                                           IFR=ifr,
                                           RMS_IFR=rms_ifr))

if __name__ == '__main__':
    print('This is not a script anymore')
