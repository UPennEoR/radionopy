import healpy as hp
import numpy as np
from matplotlib import pyplot as plt
from radiono import rm
from mpl_toolkits.basemap import Basemap
import sys,os

NSIDE=16
"""
YYYY='2011'
MM='04'
DD='11'
"""
YYYY='2012'
MM='02'
DD='13'

def IndexToDeclRa(index,deg=False):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    if deg: return -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)
    else: return theta-np.pi/2., 2.*np.pi-phi
def DeclRaToIndex(decl,RA):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(360.-RA))


IM = rm.IonoMap('30d43m17.5ss','21d25m41.9se',['%s-%s-%s'%(YYYY,MM,DD)])
#is this silly?
Plat=-30.-(43./60.)-(17.5/3600.)
Plon=21.+(25./60.)+(41.9/3600.)

ipix = range(hp.nside2npix(NSIDE))
ras,decs = [],[]
for p in ipix:
    dec,ra = IndexToDeclRa(p)
    ras.append(ra)
    decs.append(dec)

tec,rmstec,ionh = IM.ionex_data(YYYY,MM,DD)
IM.get_radec_RM(ras,decs)

mc = Basemap(projection='cyl', resolution='c')
mo = Basemap(projection='ortho',lon_0=Plon,lat_0=Plat,resolution='c')
 
for ut in range(24):
    print 'UT=%i'%ut
    
    tecUT=hp.cartview(0.1*tec[ut,:],flip='geo',max=90,min=0,return_projected_map=True)
    mc.drawcoastlines()
    _tec = mc.imshow(tecUT, alpha=0.6,vmax=90,vmin=0)
    plt.title('UT %i hr, %s-%s-%s, TEC'%(ut,YYYY,MM,DD))
    mc.colorbar(location='bottom')
    #plt.show()
    plt.savefig('TEC_%s-%s-%s_UT%i.png'%(YYYY,MM,DD,ut))
    plt.close()
    
    hp.orthview(IM.RMs[0,ut,:],half_sky=True,rot=[0,90],flip='astro',max=5,min=-2,unit=r'${\rm rad\,m^{-2}}$',return_projected_map=True)
    #mo.drawcoastlines()
    #_rm = mo.imshow(rmUT,alpha=0.6,vmax=3,vmin=-3)
    plt.title('UT %i hr, %s-%s-%s'%(ut,YYYY,MM,DD))
    #mo.colorbar(location='bottom')
    plt.savefig('RM_%s-%s-%s_UT%i.png'%(YYYY,MM,DD,ut))
    plt.close()
