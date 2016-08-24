from radiono import rm
from radiono import utils as ut
import healpy as hp, numpy as np
from matplotlib import pyplot as plt
from astropy.time import Time
# Moore et al.: 7 Dec 2011 to 27 Feb 2012
#dates = (('2011-12', range(6, 32)),('2012-01', range(1, 32)),('2012-02', range(1, 29)))

dates = []
for d in range(6,32): dates.append('-'.join(['2011-12',ut.std_hour(d)]))
for d in range(1,32): dates.append('-'.join(['2012-01',ut.std_hour(d)]))
for d in range(1,29): dates.append('-'.join(['2012-02',ut.std_hour(d)]))

# PAPER INFO
lat_str = '30d43m17.5ss'
lon_str = '21d25m41.9se'

IM = rm.IonoMap(lat_str,lon_str,dates)
ras = [np.radians(15.),np.radians(15.*4.),np.radians(15.*8.)] #1,4,8 hrs
dec = np.radians(-30.-(43./60.)-(17.5/3600.))
decs = [dec,dec,dec] #zenith

IM.get_radec_RM(ras,decs)
#import IPython;IPython.embed()

#GMST(in hours) = 6.656306 + 0.0657098242*(JD0-2445700.5) + 1.0027379093*UT
#LST = MOD [(GMST - (longitude west of Greenwich)*(24/360)),24]

ts = Time(dates,format='isot')
gmsts = np.zeros(len(dates),24)
lsts = np.zeros(len(dates,24))
for i,t in enumerate(ts):
    jd0 = ts[i].jd
    _gmsts = np.mod(6.656306 + 0.0657098242*(jd0-2445700.5) + 1.0027379093*np.linspace(0,23,num=24),24)
    gmsts[i,:] = _gmsts
    lon_wrt_W = -21.-(25./60.)-(41.9/3600.)
    lsts[i,:] = np.mod(np.mod(_gmsts - (lon_wrt_W*24./360.),24))

data = {}
for k,day in enumerate(dates)):
    data[day] = {}
    

    
