import healpy as hp
import numpy as np
from matplotlib import pyplot as plt
from radiono import rm
from mpl_toolkits.basemap import Basemap

NSIDE=16
YYYY='2011'
MM='04'
DD='11'


def IndexToDeclRa(index,deg=False):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    if deg: return -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)
    else: return theta-np.pi/2., 2.*np.pi-phi
def DeclRaToIndex(decl,RA):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(360.-RA))


IM = rm.IonoMap('0d00m0.0sn','0d00m00.0se',['%s-%s-%s'%(YYYY,MM,DD)])
#is this silly?
ipix = range(hp.nside2npix(NSIDE))
ras,decs = [],[]
for p in ipix:
    dec,ra = IndexToDeclRa(p)
    ras.append(ra)
    decs.append(dec)

tec,rmstec,ionh = IM.ionex_data(YYYY,MM,DD)
IM.get_radec_RM(ras,decs)

UT=6
tecUT=hp.cartview(tec[UT,:],flip='geo',return_projected_map=True)
rmUT = hp.cartview(IM.RMs[0,UT,:],flip='geo',return_projected_map=True)

m = Basemap(projection='cyl', resolution='c')

plt.figure(0)
m.drawcoastlines()
_tec = m.imshow(tecUT, alpha=0.6)
plt.title('UT %i hr, %s-%s-%s, TEC'%(UT,YYYY,MM,DD))
m.colorbar(location='bottom')
plt.show()
plt.close()

plt.figure(0)
m.drawcoastlines()
_rm = m.imshow(rmUT, alpha=0.6)
plt.title('UT %i hr, %s-%s-%s, RM'%(UT,YYYY,MM,DD))
m.colorbar(location='bottom')
plt.show()
plt.close()
