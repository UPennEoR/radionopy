'''
scripts.ion_altaz

purpose | script to generate RM data from IONEX file using AZs and ALTs

Functions
---------
maps2npz | writes maps to npz files
'''
from __future__ import print_function
import radiono as rad
from radiono import rm

if __name__ == '__main__':
    # PAPER INFO
    lat_str = '30d43m17.5ss'
    lon_str = '21d25m41.9se'
    
    #time_str = '2012-02-13T00:00:00'

    time_part = 'T00:00:00'
    # Moore et al.: 7 Dec 2011 to 27 Feb 2012
    """
    dates = (('2011-12', range(6, 32)),
             ('2012-01', range(1, 32)),
             ('2012-02', range(1, 29)))
    """
    # Kohn et al.: 18 Nov 2012 to 26 Mar 2013
    #dates = (('2012-11', range(18, 31)),
    #         ('2012-12', range(1, 32)),
    #         ('2013-01', range(1, 32)),
    #         ('2013-02', range(1, 29)),
    #         ('2013-03', range(1, 32)))

    dates = (('2012-11', range(18, 20)),)
    date_strs = ('-'.join((ym, rad.std_hour(day, verbose=False))) for ym, days in dates for day in days)
    time_strs = [''.join((date_str, time_part)) for date_str in date_strs]

    print(time_strs)
    RM = rm.RM(lat_str=lat_str, lon_str=lon_str, time_strs=time_strs, nside=16)
    RM.altaz()
    print(RM.RMs)
    #RM.maps_to_npz()    
