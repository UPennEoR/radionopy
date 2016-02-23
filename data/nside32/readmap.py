import numpy as np, healpy as hp
import pylab, os, sys
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

GENERATE = False
TEST = False



#2012-02-13T00:00:00 4h30m00s-30d43m17.5s CODG0440.12I.txt

karoo = EarthLocation(lat=-30.76528*u.deg,lon=21.42831*u.deg,height=1000*u.m)
datestr='2012-02-13T22:00:00'
time = Time('2012-02-13 22:00:00')
PAPERlat='30d43m17.5ss'
PAPERlon='21d25m41.9se'
IONEX='CODG0440.12I'

nside=32
npix = hp.nside2npix(nside)
ipix = np.arange(npix)
theta,phi = hp.pix2ang(nside,ipix)

alt = (90.-np.degrees(np.array(theta)))*u.degree
az = (np.degrees(np.array(phi)))*u.degree

altaz = SkyCoord(alt=alt,az=az,obstime=time,frame='altaz',location=karoo)

ra_h = altaz.icrs.ra.hms[0]
ra_m = altaz.icrs.ra.hms[1]
ra_s = altaz.icrs.ra.hms[2]

dec_d = altaz.icrs.dec.dms[0]
dec_m = altaz.icrs.dec.dms[1]
dec_s = altaz.icrs.dec.dms[2]

storage = np.zeros((npix, 15, 2)) #number of pixels, 15 UTs each, RMs and eRMs for each UT

UT0 = np.zeros((npix))
UT4 = np.zeros((npix))
UT6 = np.zeros((npix))
UT22 = np.zeros((npix))
eUT22 = np.zeros((npix))
"""
UT1_11pm
UT2_12am
UT3_1am
UT4_2am
"""
c = 0
for p in range(npix):
	
	rastr=str(int(ra_h[p]))+'h'+str(int(ra_m[p]))+'m'+str(ra_s[p])+'s'
	if dec_d[p]>0: decstr = '+'+str(int(dec_d[p]))+'d'+str(int(abs(dec_m[p])))+'m'+str(abs(dec_s[p]))+'s'
	else: decstr=str(int(dec_d[p]))+'d'+str(int(abs(dec_m[p])))+'m'+str(abs(dec_s[p]))+'s'
	
	if GENERATE:
		if TEST: 
			print 'IONFRM.py %s%s %s %s %s %s'%(rastr,decstr,PAPERlat,PAPERlon,datestr,IONEX)
		else: os.system('IONFRM.py %s%s %s %s %s %s'%(rastr,decstr,PAPERlat,PAPERlon,datestr,IONEX))
	else:
		#READ
		filename = '2012-02-13'+rastr+decstr+'IonRM.txt'
		#print 'Reading %s'%filename
		try: 
			UT, _, _, RM, eRM = np.loadtxt(filename,unpack=True)
			#print p, UT.shape
			try:
				for i, ut in enumerate(UT):
					if ut==0.: UT0[p]=RM[i]
					if ut==4.: UT4[p]=RM[i]
					if ut==6.: UT6[p]=RM[i]
					if ut==22.: 
						UT22[p]=RM[i]
						eUT22[p]=eRM[i]
			except TypeError: continue
		except IOError:
			print 'issue with %s'%filename
			UT22[p] = np.nan
			c+=1
			continue
print c,'IOErrors'	
hp.orthview(UT22,rot=[0,90],max=2,unit=r'rad m$^{-1}$',title='SAST 00:00 2012-02-13',half_sky=True)
#pylab.show()
print UT22
print 'Interpolating NaNs'
UT22_interp = np.zeros_like(UT22)
c=0
for i,val in enumerate(UT22):
	if np.isnan(val):
		c+=1
		theta,phi = hp.pix2ang(nside,i)
		neybs = hp.get_interp_weights(nside,theta,phi=phi) 
		
		v = np.nanmean(UT22[neybs[0]])
		UT22_interp[i] = v
		
	else: UT22_interp[i]=val
print c,'NaN vals'	
hp.orthview(UT22_interp,rot=[0,90],min=0,max=2,unit=r'rad m$^{-1}$',title='SAST 00:00 2012-02-13',half_sky=True)
#pylab.show()

UT22_interp_2 = np.zeros_like(UT22)
c=0
for i,val in enumerate(UT22_interp):
	if np.isnan(val):
		c+=1
		theta,phi = hp.pix2ang(nside,i)
		neybs = hp.get_interp_weights(nside,theta,phi=phi) 
		
		v = np.nanmean(UT22_interp[neybs[0]])
		
		UT22_interp_2[i] = v
		
	else: UT22_interp_2[i]=val
print c,'NaN vals'
hp.orthview(UT22_interp_2,rot=[0,90],min=0,max=2,unit=r'rad m$^{-2}$',title='SAST 00:00 2012-02-13',half_sky=True)
pylab.show()

np.savez('mapdata.npz',map=UT22_interp_2)
hp.write_map('mapdata.fits',UT22_interp_2)
