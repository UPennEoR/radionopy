from __future__ import print_function
import os
import sys
import numpy as np
import pylab as plt
#import healpy as hp
#import plotly.plotly as py
import jdcal
import rad

def std_hour(UT):
    if UT < 10:
        hour = '0{hour}'.format(hour=int(UT))
    else:
        hour = '{hour}'.format(hour=int(UT))

    return hour

def get_LST(UT, year, month, day):
    JDtup = jdcal.gcal2jd(year, month, day)
    JD = JDtup[0] + JDtup[1] + 0.5

    gmst = 6.656306 + 0.0657098242 * (JD - 2445700.5) + 1.0027379093 * UT
    lst = (gmst + 21.428 * (24. / 360.)) % 24

    return lst

base_path = os.path.expanduser('~/radionopy')
RM_dir = os.path.join(base_path, 'RM_files')

#time_strs = ('2012-02-13T00:00:00',)

time_part = 'T00:00:00'
# 7 Dec 2011 to 27 Feb 2012
dates = (('2011-12', range(6, 32)),
         ('2012-01', range(1, 32)),
         ('2012-02', range(1, 29)))
date_strs = ('-'.join((ym, std_hour(day))) for ym, days in dates for day in days)
time_strs = (''.join((date_str, time_part)) for date_str in date_strs)

#to_plot = []
lsts = ('01', '04', '08')
lst_dict = {lst: [] for lst in lsts}
for time_str in time_strs:
    year, month, day = time_str.split('T')[0].split('-')
    date = time_str.split('T')[0]

    print(time_str)
    for num in range(24):
        num = (num - 2) % 24
        lst = std_hour(get_LST(num, year, month, day))
        my_rad = os.path.join(RM_dir, '{date}/IonRM{num}.txt'.format(date=date, num=std_hour(num)))
        UT, TEC, B, RM, dRM = np.loadtxt(my_rad, unpack=True)
        RM = RM[(RM >= 0) & (RM <= 2)]
        if lst in lsts:
            lst_dict[lst].append(RM)
        #to_plot.append(RM)
        #print(lst)
#plt.show()

#plt.title("Gaussian Histogram")
#plt.xlabel("Value")
#plt.ylabel("Frequency")

for lst, RMs in lst_dict.items():
    all_RM = np.concatenate(RMs)
    #print(list(all_RM))
    plt.figure(lst)
    plt.clf()
    plt.hist(all_RM)

#for i, t in enumerate(to_plot):
#    plt.figure(i)
#    plt.clf()
#    plt.hist(t)

plt.show()
