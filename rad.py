import os
import sys
import numpy as np
import subprocess
#import pylab as plt
#import healpy as hp
from astropy import units as u
from astropy import constants as c
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle, Latitude, Longitude

# Defining some variables for further use
### Make the base path settable
base_path = os.path.expanduser('~/radionopy')
TECU = 1e16
TEC2m2 = 0.1 * TECU
earth_radius = c.R_earth.value #6371000.0 # in meters
tesla_to_gauss = 1e4

def gen_IONEX_list(IONEX_list):
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
        elif file_data.split()[-1] == 'DLON':
            start_lon, end_lon, step_lon = [float(data_item) for data_item in file_data.split()[:3]]
        elif file_data.split()[-1] == 'DLAT':
            start_lat, end_lat, step_lat = [float(data_item) for data_item in file_data.split()[:3]]

    return base_IONEX_list, RMS_IONEX_list, number_of_maps, ion_h, start_lon, end_lon, step_lon, start_lat, end_lat, step_lat

def read_IONEX_TEC(filename):
    #==========================================================================
    # Reading and storing only the TEC values of 1 day
    # (13 maps) into a 3D array

    # Opening and reading the IONEX file into memory
    with open(filename, 'r') as read_file:
        linestring = read_file.read()
        IONEX_list = linestring.split('\n')

    # creating a new array without the header and only
    # with the TEC maps
    base_IONEX_list, RMS_IONEX_list, number_of_maps, ion_h,\
    start_lon, end_lon, step_lon,\
    start_lat, end_lat, step_lat = gen_IONEX_list(IONEX_list)

    # Variables that indicate the number of points in Lat. and Lon.
    points_lon = ((end_lon - start_lon) / step_lon) + 1
    points_lat = ((end_lat - start_lat) / step_lat) + 1

    print(start_lon, end_lon, step_lon)
    print(start_lat, end_lat, step_lat)
    print(points_lon, points_lat)

    # What are the Lat/Lon coords?
    longitude = np.linspace(start_lon, end_lon, num=points_lon)
    latitude = np.linspace(start_lat, end_lat, num=points_lat)

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
    
    return TEC, RMS_TEC, (start_lon, step_lon, points_lon, start_lat, step_lat, points_lat, number_of_maps, tec_a, rms_a, ion_h * 1000.0)
    #==========================================================================

def interp(points_lat, points_lon, number_of_maps, total_maps, a):
    time_count = 1.0
    #==========================================================================================
    # producing interpolated TEC maps, and consequently a new array that will 
    # contain 25 TEC maps in total. The interpolation method used is the second
    # one indicated in the IONEX manual

    # creating a new array that will contain 25 maps in total 
    newa = np.zeros((total_maps, points_lat, points_lon))
    inc = 0
    for item in range(int(number_of_maps)):
        newa[inc, :, :] = a[item, :, :]
        inc = inc + 2

    # performing the interpolation to create 12 addional maps 
    # from the 13 TEC maps available
    time_int = int(time_count)
    while time_int <= (total_maps - 2):
        for lon in range(int(points_lon)):
            # interpolation type 2:
            # newa[int(time_count), :, lon] = 0.5 * newa[int(time_count) - 1, :, lon] + 0.5 * newa[int(time_count) + 1, :, lon]
            # interpolation type 3 ( 3 or 4 columns to the right and left of the odd maps have values of zero
            # Correct for this):
            if (lon >= 4) and (lon <= (points_lon - 4)):
                newa[time_int, :, lon] = 0.5 * newa[time_int - 1, :, lon + 3] + 0.5 * newa[time_int + 1, :, lon - 3] 
        time_int = time_int + 2

    return newa

def interp_TEC(TEC, UT, coord_lon, coord_lat, info, newa):
    start_lon, step_lon, points_lon, start_lat, step_lat, points_lat, number_of_maps, _ = info
    total_maps = 25

    #=========================================================================
    # Finding out the TEC value for the coordinates given
    # at every hour

    # Locating the 4 points in the IONEX grid map which surround
    # the coordinate you want to calculate the TEC value from  
    index_lat = 0
    index_lon = 0
    n = 0
    m = 0

    for lon in range(int(points_lon)):
        if (coord_lon > (start_lon + (n + 1) * step_lon) and coord_lon < (start_lon + (n + 2) * step_lon)):
            lower_index_lon =  n + 1
            higher_index_lon = n + 2
        n = n + 1
    for lat in range(int(points_lat)):
        if (coord_lat < (start_lat + (m + 1) * step_lat) and coord_lat > (start_lat + (m + 2) * step_lat)):
            lower_index_lat =  m + 1
            higher_index_lat = m + 2
        m = m + 1

    # Using the 4-point formula indicated in the IONEX manual
    # The TEC value at the coordinates you desire for every 
    # hour are estimated 
    diff_lon = coord_lon - (start_lon + lower_index_lon * step_lon)
    p = diff_lon / step_lon
    diff_lat = coord_lat - (start_lat + lower_index_lat * step_lat)
    q = diff_lat / step_lat
    TEC_values = []
    for m in range(total_maps):
        TEC_values.append((1.0 - p) * (1.0 - q) * newa[m, lower_index_lat, lower_index_lon]\
                          + p * (1.0 - q) * newa[m, lower_index_lat, higher_index_lon]\
                          + q * (1.0 - p) * newa[m, higher_index_lat, lower_index_lon]\
                          + p * q * newa[m, higher_index_lat, higher_index_lon])
    #=========================================================================

    #return {'TEC_values': np.array(TEC_values), 'a': np.array(a), 'newa': np.array(newa)}
    return np.array(TEC_values)[UT]

def punct_ion_offset(lat_obs, az_src, zen_src, ion_height):
    #earth_radius = 6371000.0 # in meters

    # The 2-D sine rule gives the zenith angle at the
    # Ionospheric piercing point
    zen_punct = np.arcsin((earth_radius * np.sin(zen_src)) / (earth_radius + ion_height)) 

    # Use the sum of the internal angles of a triange to determine theta
    theta = zen_src - zen_punct

    # The cosine rule for spherical triangles gives us the latitude
    # at the IPP
    lat_ion = np.arcsin(np.sin(lat_obs) * np.cos(theta) + np.cos(lat_obs) * np.sin(theta) * np.cos(az_src)) 
    off_lat = lat_ion - lat_obs # latitude difference

    # Longitude difference using the 3-D sine rule (or for spherical triangles)
    off_lon = np.arcsin(np.sin(az_src) * np.sin(theta) / np.cos(lat_ion))

    # Azimuth at the IPP using the 3-D sine rule
    s_az_ion = np.sin(az_src) * np.cos(lat_obs) / np.cos(lat_ion)
    az_punct = np.arcsin(s_az_ion)

    return off_lon, off_lat, az_punct, zen_punct

def get_coords(lon_str, lat_str, lon_obs, lat_obs, off_lon, off_lat):
    if lon_str[-1] == 'e':
        lon_val = 1
    elif lon_str[-1] == 'w':
        lon_val = -1
    if lat_str[-1] == 's':
        lat_val = -1
    elif lat_str[-1] == 'n':
        lat_val = 1

    coord_lon = lon_val * (lon_obs.value + off_lon)
    coord_lat = lat_val * (lat_obs.value + off_lat)

    return coord_lon, coord_lat

def TEC_paths(TEC, RMS_TEC, UT, coord_lon, coord_lat, zen_punct, info, rms_info, newa, rmsa):
    VTEC = interp_TEC(TEC, UT, coord_lon, coord_lat, info, newa)
    TEC_path = VTEC * TEC2m2 / np.cos(zen_punct) # from vertical TEC to line of sight TEC

    VRMS_TEC = interp_TEC(RMS_TEC, UT, coord_lon, coord_lat, rms_info, rmsa)
    RMS_TEC_path = VRMS_TEC * TEC2m2 / np.cos(zen_punct) # from vertical RMS_TEC to line of sight RMS_TEC

    return TEC_path, RMS_TEC_path

def B_IGRF(year, month, day, coord_lon, coord_lat, ion_height, az_punct, zen_punct):
    # Calculation of TEC path value for the indicated 'hour' and therefore 
    # at the IPP

    input_file = os.path.join(base_path, 'IGRF/geomag70_linux/input.txt')
    output_file = os.path.join(base_path, 'IGRF/geomag70_linux/output.txt')

    #uses lat_val, lon_val from above
    # Calculation of the total magnetic field along the line of sight at the IPP
    with open(input_file, 'w') as f:
        f.write('{year},{month},{day} C K{sky_rad} {ipp_lat} {ipp_lon}'.format(year=year, month=month, day=day,
                                                                               sky_rad=(earth_radius + ion_height) / 1000.0,
                                                                               ipp_lon=coord_lon, ipp_lat=coord_lat))

    #XXX runs the geomag exe script
    script_name = os.path.join('./', base_path, 'IGRF/geomag70_linux/geomag70')
    script_data = os.path.join(base_path, 'IGRF/geomag70_linux/IGRF11.COF')
    script_option = 'f'
    subprocess.call([script_name, script_data, script_option, input_file, output_file])

    with open(output_file, 'r') as g:
        data = g.readlines()


        x_field, y_field, z_field = [abs(float(field_data)) * 1e-9 * tesla_to_gauss for field_data in data[1].split()[10:13]]
        tot_field = z_field * np.cos(zen_punct) +\
                    y_field * np.sin(zen_punct) * np.sin(az_punct) +\
                    x_field * np.sin(zen_punct) * np.cos(az_punct)

    #remove files once used
    #os.remove(input_file)
    #os.remove(output_file)

    return tot_field

def get_results(lat_obs, lon_obs, alt_src, az_src, zen_src, ion_height, TEC, RMS_TEC, info, rms_info, UT, newa, rmsa):
    if (alt_src.degree > 0):
        print(i, alt_src, az_src)
        # Calculate the ionospheric piercing point.  Inputs and outputs in radians
        off_lon, off_lat, az_punct, zen_punct = punct_ion_offset(lat_obs.radian, az_src.radian, zen_src.to(u.radian).value, ion_height)
        print(off_lon, off_lat, az_punct, zen_punct)

        coord_lon, coord_lat = get_coords(lon_str, lat_str, lon_obs, lat_obs, off_lon * 180 / np.pi, off_lat * 180 / np.pi)

        TEC_path, RMS_TEC_path = TEC_paths(TEC, RMS_TEC, UT, coord_lon, coord_lat, zen_punct, info, rms_info, newa, rmsa)
        tot_field = B_IGRF(year, month, day, coord_lon, coord_lat, ion_height, az_punct, zen_punct)

        # Saving the Ionosheric RM and its corresponding
        # rms value to a file for the given 'hour' value
        IFR = 2.6e-17 * tot_field * TEC_path
        RMS_IFR = 2.6e-17 * tot_field * RMS_TEC_path

        with open(os.path.join(base_path, 'IonRM.txt'), 'a') as f:
            f.write('{hour} {TEC_path} {tot_field} {IFR} {RMS_IFR}\n'.format(hour=hour, TEC_path=TEC_path, tot_field=tot_field,
                                                                             IFR=IFR, RMS_IFR=RMS_IFR))

        return {'TEC': TEC, 'RMS_TEC': RMS_TEC, 'RM': IFR, 'RMS_RM': RMS_IFR, 'tot_field': tot_field}

if __name__ == '__main__':
    with open(os.path.join(base_path, 'IonRM.txt'), 'w') as f:
        pass
### Here we need to accept an array of RA/Dec which correspond to the
### centers of healpix pixels, and all subsequent operations should
### allow ra/dec to be vectorized

### The returned value of the function is then an RA/Dec map of the RM
### above the array

### npix = hp.nside2npix(nside)
### ipix = np.arange(npix)
### ra, dec = hp.pix2ang(nside, ipix)
### ra_dec = SkyCoord(ra=ra, dec=dec) # go from healpix theta, phi radians to astropy ra, dec
### alt_src = ra_dec.altaz.alt
### az_src = ra_dec.altaz.az
### that passes in to the new function

    ## Nominally try to reproduce the output of this command
    ## ionFRM.py 16h50m04.0s+79d11m25.0s 52d54m54.64sn 6d36m16.04se 2004-05-19T00:00:00 CODG1400.04I
    ## Echo back what he has ... 
    ra_str = '16h50m04.0s'
    dec_str = '+79d11m25.0s'
    lon_str = '6d36m16.04se'
    lat_str = '52d54m54.64sn'
    time_str = '2004-05-19T00:00:00' # This will actually work as input to the astropy Time function
    IONEX_file = 'CODG1400.04I'
    IONEX_name = os.path.join(base_path, IONEX_file)

    year, month, day = time_str.split('T')[0].split('-')

    lon_obs = Longitude(Angle(lon_str[:-1]))
    lat_obs = Latitude(Angle(lat_str[:-1]))

    location = EarthLocation(lon=lon_obs, lat=lat_obs, height=0 * u.m)
    start_time = Time(time_str)

    # Create a sky coordinate object, from which we can subsequently derive the necessary alt/az
    ra_dec = SkyCoord(ra=ra_str, dec=dec_str, location=location, obstime=start_time)

    TEC, RMS_TEC, all_info = read_IONEX_TEC(IONEX_name)

    info = all_info[:7] + (all_info[7],)
    rms_info = all_info[:7] + (all_info[8],)
    ion_height = all_info[9]

    _, _, points_lon, _, _, points_lat, number_of_maps, a = info
    _, _, _, _, _, _, _, rms_a = rms_info

    newa = interp(points_lat, points_lon, number_of_maps, 25, a)
    rmsa = interp(points_lat, points_lon, number_of_maps, 25, rms_a)

    # predict the ionospheric RM for every hour within a day 
    UTs = np.linspace(0, 23, num=24)
    
    results_dict = {}
    
    for i, UT in enumerate(UTs):
        print(UT)
        if UT < 10:
            hour = '0{hour}'.format(hour=int(UT))
        else:
            hour = '{hour}'.format(hour=int(UT))
        
        ra_dec = SkyCoord(ra=ra_str, dec=dec_str, location=location, obstime=start_time + UT * u.hr)

        # Calculate alt and az
        alt_src = ra_dec.altaz.alt
        az_src = ra_dec.altaz.az

        # zen_src is a different kind of object than Alt/Az
        zen_src = ra_dec.altaz.zen

        #thing = {'TEC': TEC, 'RMS_TEC': RMS_TEC, 'IFR': IFR, 'RMS_IFR': RMS_IFR, 'tot_field': tot_field}
        thing = get_results(lat_obs, lon_obs, alt_src, az_src, zen_src, ion_height, TEC, RMS_TEC, info, rms_info, UT, newa, rmsa)
