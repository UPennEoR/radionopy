'''
scripts.hist

purpose | script to plot RM data by LST over range of dates into histogram

Functions
---------
get_LST | generates LST from date
get_RM | finds RM data for particular date and LST
'''
from __future__ import print_function
import os
import sys
import numpy as np
import pylab as plt
from scipy import spatial
import jdcal
import radiono as rad

def get_LST(UT, year, month, day):
    '''
    gets LST for particular date and time

    Parameters
    ----------
    UT | int: numbered hour of time
    year | int: year
    month | int: numbered month of the year
    day | int: numbered day of the month

    Returns
    -------
    int: local sidereal time by hour
    '''
    JDtup = jdcal.gcal2jd(year, month, day)
    JD = JDtup[0] + JDtup[1] + 0.5

    gmst = 6.656306 + 0.0657098242 * (JD - 2445700.5) + 1.0027379093 * UT
    lst = (gmst + 21.428 * (24. / 360.)) % 24

    return lst

def get_RM(num, date, LST):
    '''
    gets RM data for particular date and LST

    Parameters
    ----------
    num | int: integer value for hour
    date | str: string representation of date
    LST | int: local sidereal time by hour
    '''
    rm_file = os.path.join(RM_dir, '{date}/IonRM{num}.txt'.format(date=date, num=rad.std_hour(num)))
    _, _, _, RM, dRM = np.loadtxt(rm_file, unpack=True)

    ra = LST * 15

    radec_file = os.path.join(RM_dir, '{date}/radec{num}.txt'.format(date=date, num=rad.std_hour(num)))
    RADECs = np.loadtxt(radec_file)
    distance, idx = spatial.KDTree(RADECs).query((ra,dec))
    #print(distance, idx, RADECs[idx])
    return RM[idx]

if __name__ == '__main__':
    dec = -30.76528

    RM_dir = os.path.join(rad.base_path, 'RM_files')

    time_part = 'T00:00:00'
    # 7 Dec 2011 to 27 Feb 2012
    dates = (('2011-12', range(6, 32)),
             ('2012-01', range(1, 32)),
             ('2012-02', range(1, 29)))
    date_strs = ('-'.join((ym, rad.std_hour(day))) for ym, days in dates for day in days)
    time_strs = (''.join((date_str, time_part)) for date_str in date_strs)

    lsts = ('01', '04', '08')
    lst_dict = {lst: [] for lst in lsts}
    #lst_dict = {rad.std_hour(lst): [] for lst in range(24)}
    for time_str in time_strs:
        year, month, day = time_str.split('T')[0].split('-')
        date = time_str.split('T')[0]

        print(time_str)
        for num in range(24):
            num = (num - 2) % 24
            LST = get_LST(num, year, month, day)
            lst = rad.std_hour(LST)
            RM = get_RM(num, date, LST)
            #RM = RM[(RM >= 0) & (RM <= 2)]
            if not 0 <= RM <= 2:
                continue
            if lst in lsts:
                lst_dict[lst].append(RM)

    for lst, RMs in lst_dict.items():
        #all_RM = np.concatenate(RMs)
        all_RM = RMs
        #print(list(all_RM))
        plt.figure(lst)
        plt.clf()
        plt.hist(all_RM, bins=20, range=[0,2])

    plt.show()
