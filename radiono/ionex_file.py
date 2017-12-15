'''
radiono.ionex_file

purpose | Module used to gather information from IONEX files

Functions
---------
pull_IONEX_file | finds correct IONEX file and uncompresses it if necessary
ftp_IONEX_file | downloads IONEX file from ftp server if not on local machine
parse_IONEX_file | pulls specific info from IONEX file into a tuple of various data
get_IONEX_data | Generate information from IONEX file to be used in further calculations
'''
from __future__ import print_function
import os
import datetime
import ftplib
import subprocess
import numpy as np
import radiono as rad

def pull_IONEX_file(year, month, day, ionex_dir=rad.ionex_dir, verbose=False):
    '''
    pulls correct IONEX file for date input
    decompresses if necessary

    Parameters
    ----------
    year | int: year
    month | int: numbered month of the year
    day | int: numbered day of the month
    ionex_dir | Optional[str]: directory in which ionex files are / will be located

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
    if verbose: print(ionex_file,ionex_file_z)

    if not os.path.exists(os.path.join(ionex_dir, ionex_file)): # or not os.path.exists(os.path.join(ionex_dir, ionex_file_z)):
        if verbose: print('Getting ionex file %s'%ionex_file_z)
        ionex_file_z = ftp_IONEX_file(year, month, day, ionex_file)
        if verbose: print('Uncompressing ionex file %s -> %s'%(ionex_file_z,ionex_file))
        os.system('gunzip --quiet '+ionex_file_z)
        if verbose: print('Done')
    return os.path.join(ionex_dir, ionex_file)

def ftp_IONEX_file(year, month, day, IONEX_file, ionex_dir=rad.ionex_dir):
    '''
    downloads IONEX file from ftp server

    Parameters
    ----------
    year | int: year
    month | int: numbered month of the year
    day | int: numbered day of the month
    IONEX_file | str: name of IONEX file to be downloaded
    ionex_dir | Optional[str]: directory in which ionex files are / will be located

    Returns
    -------
    str: name of IONEX file pulled from ftp server

    '''
    server = 'ftp.aiub.unibe.ch'
    year = str(year)
    ftp_dir = os.path.join('CODE/', year)
    IONEX_file_Z = ''.join((os.path.basename(IONEX_file), '.Z'))

    if not os.path.exists(ionex_dir):
        os.mkdir(ionex_dir)

    IONEX_file_X = os.path.join(rad.ionex_dir, os.path.basename(IONEX_file_Z))

    getting_file_str = 'Retrieving {IONEX_file_Z} for {day} {month} {year}'.format(IONEX_file_Z=IONEX_file_Z, day=day, month=month, year=year)
    print(getting_file_str)

    try:
        ftp = ftplib.FTP(server, 'anonymous', 'jaguirre@sas.upenn.edu')
        ftp.cwd(ftp_dir)
        ftp_cmd=' '.join(('RETR', IONEX_file_Z))
        ftp.retrbinary(ftp_cmd, open(IONEX_file_X, 'wb').write)
        ftp.quit()
    except:
        print('No file exists?')
        os.remove(IONEX_file_X)

    return IONEX_file_X

def parse_IONEX_file(IONEX_file):
    '''
    parse data from IONEX file

    Parameters
    ----------
    IONEX_file | str: IONEX filename

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
    # Opening and reading the IONEX file into memory
    with open(IONEX_file, 'r') as read_file:
        linestring = read_file.read()
        IONEX_list = linestring.split('\n') #this is the IONEX file in a python list

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
            number_of_maps = int(float(file_data.split()[0]))

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

def get_IONEX_data(filename, verbose=False):
    '''

    Generate information from IONEX file to be used in further calculations

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

    # IONEX_list = pull_IONEX_file(filename)

    # creating a new array without the header and only
    # with the TEC maps
    base_IONEX_list, RMS_IONEX_list, number_of_maps, ion_h,\
    start_lat, end_lat, step_lat,\
    start_lon, end_lon, step_lon = parse_IONEX_file(filename)

    # Variables that indicate the number of points in Lat. and Lon.
    points_lat = int( ((end_lat - start_lat) / step_lat) + 1 )
    points_lon = int( ((end_lon - start_lon) / step_lon) + 1 )
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
