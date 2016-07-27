'''
scripts.ion_radec

purpose | script to generate RM data from IONEX file using RAs and DECs
'''
from __future__ import print_function
import numpy as np
from astropy.coordinates import SkyCoord
from radiono import rm

if __name__ == '__main__':
    ra_strs = ('16h50m04.0s',)
    dec_strs = ('+79d11m25.0s',)
    lat_str = '52d54m54.64sn'
    lon_str = '6d36m16.04se'
    time_strs = ('2004-05-19T00:00:00',)

    radec = SkyCoord(ra=np.array(ra_strs), dec=np.array(dec_strs))
    ra, dec = radec.ra, radec.dec
    RM = rm.RM(lat_str, lon_str, time_strs)
    RM.radec(ra, dec)
    print(RM.RMs)
    print(RM.B_para)