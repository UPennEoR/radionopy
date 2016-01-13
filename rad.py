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

base_path = os.path.expanduser('~/radionopy')

def read_IONEX_TEC(filename, rms=False):
    #==========================================================================
    # Reading and storing only the TEC values of 1 day
    # (13 maps) into a 3D array

    # Opening and reading the IONEX file into memory
    with open(filename, 'r') as read_file:
        linestring = read_file.read()
        IONEX_list = linestring.split('\n')

    # creating a new array without the header and only
    # with the TEC maps
    add = 0 
    new_IONEX_list = []
    for file_data in IONEX_list[:-1]:
        if not file_data:
            continue
        if file_data.split()[-2:] == ['RMS', 'MAP']:
            add = 0
            if rms:
                add = 1
        elif file_data.split()[-2:] == ['IN', 'FILE']:
            number_of_maps = float(file_data.split()[0])

        if not rms:
            if file_data.split()[0] == 'END' and file_data.split()[2] == 'HEADER':
                add = 1

        if add == 1:
            new_IONEX_list.append(file_data)

        if file_data.split()[-1] == 'DHGT':
            ion_h = float(file_data.split()[0])
        elif file_data.split()[-1] == 'DLON':
            start_lon, end_lon, step_lon = [float(data_item) for data_item in file_data.split()[:3]]
        elif file_data.split()[-1] == 'DLAT':
            start_lat, end_lat, step_lat = [float(data_item) for data_item in file_data.split()[:3]]

    # Variables that indicate the number of points in Lat. and Lon.
    points_lon = ((end_lon - start_lon) / step_lon) + 1
    points_lat = ((end_lat - start_lat) / step_lat) + 1

    print(start_lon, end_lon, step_lon)
    print(start_lat, end_lat, step_lat)
    print(points_lon, points_lat)

    # What are the Lat/Lon coords?
    longitude = np.linspace(start_lon, end_lon, num=points_lon)
    latitude = np.linspace(start_lat, end_lat, num=points_lat)

    # 3D array that will contain TEC values only
    a = np.zeros((number_of_maps, points_lat, points_lon))

    # Selecting only the TEC values to store in the 3-D array
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

    TEC =  {'TEC': np.array(a), 'lat': latitude, 'lon': longitude, 'AltIon': ion_h * 1000.0}
    return TEC, (start_lon, step_lon, points_lon, start_lat, step_lat, points_lat, number_of_maps, a)
    #==========================================================================

def interp_TEC(TEC, UT, coord_lon, coord_lat, info):
    total_maps = 25

    start_lon, step_lon, points_lon, start_lat, step_lat, points_lat, number_of_maps, a = info

    time_count = 1.0
    #new_UT = UT
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
    while int(time_count) <= (total_maps - 2):
        for lat in range(int(points_lat)):
            for lon in range(int(points_lon)):
                # interpolation type 2:
                # newa[int(time_count), lat, lon] = 0.5 * newa[int(time_count) - 1, lat, lon] + 0.5 * newa[int(time_count) + 1, lat, lon]
                # interpolation type 3 ( 3 or 4 columns to the right and left of the odd maps have values of zero
                # Correct for this):
                if (lon >= 4) and (lon <= (points_lon - 4)):
                    newa[int(time_count), lat, lon] = 0.5 * newa[int(time_count) - 1, lat, lon + 3] + 0.5 * newa[int(time_count) + 1, lat, lon - 3] 
        time_count = time_count + 2.0
    #==========================================================================================


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

def punct_ion_offset(lat_obs, az_source, zen_source, alt_ion):
    #earth_radius = 6371000.0 # in meters

    # The 2-D sine rule gives the zenith angle at the
    # Ionospheric piercing point
    zen_punct = np.arcsin((earth_radius * np.sin(zen_source)) / (earth_radius + alt_ion)) 

    # Use the sum of the internal angles of a triange to determine theta
    theta = zen_source - zen_punct

    # The cosine rule for spherical triangles gives us the latitude
    # at the IPP
    lat_ion = np.arcsin(np.sin(lat_obs) * np.cos(theta) + np.cos(lat_obs) * np.sin(theta) * np.cos(az_source)) 
    d_lat = lat_ion - lat_obs # latitude difference

    # Longitude difference using the 3-D sine rule (or for spherical triangles)
    d_lon = np.arcsin(np.sin(az_source) * np.sin(theta) / np.cos(lat_ion))

    # Azimuth at the IPP using the 3-D sine rule
    s_az_ion = np.sin(az_source) * np.cos(lat_obs) / np.cos(lat_ion)
    az_punct = np.arcsin(s_az_ion)

    return d_lon, d_lat, az_punct, zen_punct

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

def TEC_paths(TEC, UT, coord_lon, coord_lat, zen_punct, info):
    VTEC = interp_TEC(TEC, UT, coord_lon, coord_lat, info)
    TEC_path = VTEC * TEC2m2 / np.cos(zen_punct) # from vertical TEC to line of sight TEC

    return TEC_path

def B_IGRF(year, month, day, coord_lon, coord_lat, alt_ion, az_punct, zen_punct):
    # Calculation of TEC path value for the indicated 'hour' and therefore 
    # at the IPP

    input_file = os.path.join(base_path, 'IGRF/geomag70_linux/input.txt')
    output_file = os.path.join(base_path, 'IGRF/geomag70_linux/output.txt')

    #uses lat_val, lon_val from above
    # Calculation of the total magnetic field along the line of sight at the IPP
    with open(input_file, 'w') as f:
        f.write('{year},{month},{day} C K{sky_rad} {ipp_lat} {ipp_lon}'.format(year=year, month=month, day=day,
                                                                               sky_rad=(earth_radius + alt_ion) / 1000.0,
                                                                               ipp_lon=coord_lon, ipp_lat=coord_lat))

    #XXX runs the geomag exe script
    script_name = os.path.join('./', base_path, 'IGRF/geomag70_linux/geomag70')
    script_data = os.path.join(base_path, 'IGRF/geomag70_linux/IGRF11.COF')
    script_option = 'f'
    subprocess.call([script_name, script_data, script_option, input_file, output_file])

    with open(output_file, 'r') as g:
        data = g.readlines()


        x_field, y_field, z_field = [abs(float(field_data)) * pow(10, -9) * tesla_to_gauss for field_data in data[1].split()[10:13]]
        tot_field = z_field * np.cos(zen_punct) +\
                    y_field * np.sin(zen_punct) * np.sin(az_punct) +\
                    x_field * np.sin(zen_punct) * np.cos(az_punct)

    #remove files once used
    #os.remove(input_file)
    #os.remove(output_file)

    return tot_field

if __name__ == '__main__':
    with open(os.path.join(base_path, 'IonRM.txt'), 'w') as f:
        pass
    # Defining some variables for further use
    TECU = pow(10, 16)
    TEC2m2 = 0.1 * TECU
    earth_radius = c.R_earth.value #6371000.0 # in meters
    tesla_to_gauss = pow(10, 4)

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

    TEC, info = read_IONEX_TEC(IONEX_name)
    RMS_TEC, rms_info = read_IONEX_TEC(IONEX_name, rms=True)

    # Reading the altitude of the Ionosphere in km (from IONEX file)
    alt_ion = TEC['AltIon']

    # predict the ionospheric RM for every hour within a day 
    UTs = np.linspace(0, 23, num=24)
    for i, UT in enumerate(UTs):
        print(UT)
        if UT < 10:
            hour = '0{hour}'.format(hour=UT)
        else:
            hour = '{hour}'.format(hour=UT)
        
        ra_dec = SkyCoord(ra=ra_str, dec=dec_str, location=location, obstime=start_time + UT * u.hr)

        # Calculate alt and az
        alt_source = ra_dec.altaz.alt
        az_source = ra_dec.altaz.az

        # zen_source is a different kind of object than Alt/Az
        zen_source = ra_dec.altaz.zen

        if (alt_source.degree > 0):
            print(i, alt_source, az_source)
            # Calculate the ionospheric piercing point.  Inputs and outputs in radians
            off_lon, off_lat, az_punct, zen_punct = punct_ion_offset(lat_obs.radian, az_source.radian, zen_source.to(u.radian).value, alt_ion)
            print(off_lon, off_lat, az_punct, zen_punct)

            coord_lon, coord_lat = get_coords(lon_str, lat_str, lon_obs, lat_obs, off_lon, off_lat)

            TEC_path = TEC_paths(TEC, UT, coord_lon, coord_lat, zen_punct, info)
            RMS_TEC_path = TEC_paths(RMS_TEC, UT, coord_lon, coord_lat, zen_punct, rms_info)
            tot_field = B_IGRF(year, month, day, coord_lon, coord_lat, alt_ion, az_punct, zen_punct)

            # Saving the Ionosheric RM and its corresponding
            # rms value to a file for the given 'hour' value
            IFR = 2.6 * pow(10, -17) * tot_field * TEC_path
            RMS_IFR = 2.6 * pow(10, -17) * tot_field * RMS_TEC_path
            with open(os.path.join(base_path, 'IonRM.txt'), 'a') as f:
                f.write('{hour} {TEC_path} {tot_field} {IFR} {RMS_IFR}\n'.format(hour=hour, TEC_path=TEC_path, tot_field=tot_field,
                                                                                 IFR=IFR, RMS_IFR=RMS_IFR))
