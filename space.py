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
