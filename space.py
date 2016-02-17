def find_nearest(point, vector):
    diff = np.abs(vector - point)
    wh = np.where(diff == diff.min())
    return wh[0]

def interp_space(TEC, UT, coord_lat, coord_lon, info, newa):
    start_lat, step_lat, points_lat, start_lon, step_lon, points_lon, number_of_maps, _ = info
    total_maps = 25

    latitudes = TEC['lat']
    longitudes = TEC['lon']

    #=========================================================================
    # Finding out the TEC value for the coordinates given
    # at every hour

    # Locating the 4 points in the IONEX grid map which surround
    # the coordinate you want to calculate the TEC value from  
    m = 0
    n = 0

    lower_index_lat = (find_nearest(coord_lat, latitudes) + 1) % len(latitudes)
    lower_index_lon = (find_nearest(coord_lon, longitudes) + 1) % len(longitudes)
    higher_index_lat = (lower_index_lat + 1) % len(latitudes)
    higher_index_lon = (lower_index_lon + 1) % len(longitudes)

    # Using the 4-point formula indicated in the IONEX manual
    # The TEC value at the coordinates you desire for every 
    # hour are estimated 
    diff_lat = coord_lat - (start_lat + lower_index_lat * step_lat)
    q = diff_lat / step_lat
    diff_lon = coord_lon - (start_lon + lower_index_lon * step_lon)
    p = diff_lon / step_lon
    TEC_values = []
    try:
        for m in range(total_maps):
            TEC_values.append(
                (1.0 - p) * (1.0 - q) * newa[m, lower_index_lat, lower_index_lon]\
                + p * (1.0 - q) * newa[m, lower_index_lat, higher_index_lon]\
                + q * (1.0 - p) * newa[m, higher_index_lat, lower_index_lon]\
                + p * q * newa[m, higher_index_lat, higher_index_lon])
    except:
        print(lower_index_lat)
        print(lower_index_lon)
        print(higher_index_lat)
        print(higher_index_lon)
    #=========================================================================

    return np.array(TEC_values)[UT][0]

    #newa = interp_time(points_lat, points_lon, number_of_maps, 25, a)
    #rmsa = interp_time(points_lat, points_lon, number_of_maps, 25, rms_a)
def interp_time(points_lat, points_lon, number_of_maps, total_maps, a):
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
            #if (lon >= 4) and (lon <= (points_lon - 4)):
            #    newa[time_int, :, lon] = 0.5 * newa[time_int - 1, :, lon + 3] + 0.5 * newa[time_int + 1, :, lon - 3] 
            newa[time_int, :, lon] = 0.5 * newa[time_int - 1, :, (lon + 3) % int(points_lon)] + 0.5 * newa[time_int + 1, :, lon - 3] 
        time_int = time_int + 2

    return newa
