import os, sys
import numpy as np
from radiono import rm
import matplotlib.pyplot as plt

RA = '23h23m27.9s'
Dec='+58d48m42.4s'
lat='53d00m00.00sn'
lng='6d00m00.00se' 

year='2011'
month='04'
day='11'

ionexfile = 'CODG1010.11I'

time = year+'-'+month+'-'+day+'T00:00:00'

ionFRMcommand = 'ionFRM.py '+RA+Dec+' '+lat+' '+lng+' '+time+' '+'IONEX/Data/'+ionexfile
ionfrfile = '2011-04-1123h23m27.9s+58d48m42.4sIonRM.txt'
if not os.path.exists(ionfrfile): os.system(ionFRMcommand)
ut,tec,absB,RM,dRM = np.loadtxt(ionfrfile,unpack=True)

RA_d = 23.+(23./60.)+(27.9/3600.)
RA_d *= 15.
Dec_d = 58.+(48./60.)+(42.4/3600.)

IM = rm.IonoMap(lat,lng,['%s-%s-%s'%(year,month,day)])
IM.calc_radec_rm([np.radians(RA_d)],[np.radians(Dec_d)])

plt.errorbar(ut,RM,yerr=dRM,fmt='ro',label=r'${\tt ionFR}$')
plt.errorbar(range(24),IM.RMs[0,:],yerr=IM.dRMs[0,:],fmt='bo',label=r'${\tt radionopy}$')
plt.suptitle('2011-04-11')
plt.xlabel('UT [hours]')
plt.ylabel(r'$\rm\phi_{iono}\,[rad\,m^{-2}]$')
plt.legend(loc=4)
plt.xlim(0,23)
plt.ylim(0,2.5)
plt.show()