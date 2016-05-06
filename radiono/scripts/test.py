from __future__ import print_function
import numpy as np
import radiono as rad
from radiono import rm

ra_strs = ('16h50m04.0s',)
dec_strs = ('+79d11m25.0s',)
lat_str = '52d54m54.64sn'
lon_str = '6d36m16.04se'
time_strs = ('2004-05-19T00:00:00',)

UTs = np.linspace(0, 23, num=24)

test = rm.RM(lat_str, lon_str, time_strs)
#test._radec(ra_strs, dec_strs, UTs)
test._altaz()
print(test.RMs)
