from radiono import rm
from radiono import utils as ut
import healpy as hp, numpy as np, os
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

if not os.path.exists('./1h4h8h.npz'):
    print 'getting RMs'
    IM.get_radec_RM(ras,decs)
    np.savez('./1h4h8h.npz',RM=IM.RMs,dRM=IM.dRMs,dates=dates)
else:
    d = np.load('./1h4h8h.npz')
    dvec = np.zeros( (len(d['dates']),3,2) )
    datevec = []
    for i,ra in enumerate([15.,60.,120.]):#i counts RA pointings, in degrees
        for j,dt in enumerate(d['dates']):#j counts date in season
            nT = ut.nextTransit(ut.ion2ephDate(dt),ra,dec)
            pTB = ut.parseTransitBasic(nT,SunCheck=True)
            
            #its next transit the next day?
            
            if not map(int,dt.split('-'))==map(int,pTB[0].split('-')) and pTB[1]<24:
                J=j+1
                k = pTB[1]
            elif pTB[2] and pTB[1]<24: 
                zeroRM=True #is the sun up?
                J = j
                k = pTB[1]
            elif pTB[1]>=24:
                J = j+1
                k = pTB[1]%24
            else:
                zeroRM=False
                J = j 
                k = pTB[1]
            rm,drm = d['RM'][J,k,i],d['dRM'][J,k,i]
            if not zeroRM: dvec[J,i,:] = [rm,drm]
            else: dvec[J,i,:] = [np.nan,np.nan]
            datevec.append([ra,dt,nT])

f,axarr = plt.subplots(3,1,sharex=True)
for N in range(3):
    dhist = dvec[:,N,0][~np.isnan(dvec[:,N,0])]
    bw=0.3 #bin width in rad/m^2
    axarr[N].hist(dhist,bins=np.arange(min(dhist),max(dhist)+bw,bw))
    axarr[N].set_ylim(0,60)
axarr[2].set_xlim(-2,7)
plt.show()

import IPython;IPython.embed()

