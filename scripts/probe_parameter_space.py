import healpy as hp, numpy as np, matplotlib.pyplot as plt, sys, os, aipy, pickle
from bm_prms import prms
from radiono import rm
from radiono import utilts as ut

NSIDE=16

##### These should be moved into radiono/something.py (utils or rm)
def IndexToDeclRa(index,deg=False):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    if deg: return -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)
    else: return theta-np.pi/2., 2.*np.pi-phi
def DeclRaToIndex(decl,RA):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(360.-RA))
def beamPAPER(freq_MHz,nside=512,lmax=20,mmax=20):
    """
    Construct a normalized AIPY beam for a single frequency (XX-pol)
    """
    freq=freq_MHz/1000. #aipy wants GHz
    bm = prms['beam'](np.array([freq]),nside=nside,lmax=lmax,mmax=mmax,deg=7)
    bm.set_params(prms['bm_prms'])
    px = range(hp.nside2npix(nside))
    xyz = hp.pix2vec(nside,px)
    poly = np.array([h.map[px] for h in bm.hmap])
    Axx = np.polyval(poly,freq)
    Axx = np.where(xyz[-1] >= 0, Axx, 0)
    Axx /= Axx.max()
    Axx = Axx*Axx
    beam = rotate_hmap(Axx,[0,-90])
    return beam

def rotate_hmap(map,rot):
	npix = map.shape[0]
	nside = hp.npix2nside(npix)

	rotmap = np.zeros(npix)
	ipix = np.arange(npix)
	t,p = hp.pix2ang(nside,ipix)

	r = hp.Rotator(rot=rot)

	# For each pixel in the new map, find where it would have come 
	# from in the old    
	trot,prot = r(t,p)
	ipix_rot = hp.ang2pix(nside,trot,prot)

	rotmap = map[ipix_rot]

	return rotmap
	
Plat=-30.-(43./60.)-(17.5/3600.)
Plon=21.+(25./60.)+(41.9/3600.)

ipix = range(hp.nside2npix(NSIDE))
ras,decs = [],[]
for p in ipix:
    dec,ra = IndexToDeclRa(p)
    ras.append(ra)
    decs.append(dec)
###### 

#dates = ['2012-02-10', '2012-02-11', '2012-02-13', '2012-02-14', '2012-02-15', '2012-02-16', '2012-02-17','2012-02-18', '2012-02-19', '2012-02-20', '2012-02-21',  '2012-02-22']
dates = ['2011-12-06', '2011-12-07', '2011-12-08', '2011-12-09', '2011-12-10', '2011-12-11', '2011-12-12', '2011-12-13', '2011-12-14', '2011-12-15', '2011-12-16', '2011-12-17', '2011-12-18', '2011-12-19', '2011-12-20', '2011-12-21', '2011-12-22', '2011-12-23', '2011-12-24', '2011-12-25', '2011-12-26', '2011-12-27', '2011-12-28', '2011-12-29', '2011-12-30', '2011-12-31', '2012-01-01', '2012-01-02', '2012-01-03', '2012-01-04', '2012-01-05', '2012-01-06', '2012-01-07', '2012-01-08', '2012-01-09', '2012-01-10', '2012-01-11', '2012-01-12', '2012-01-13', '2012-01-14', '2012-01-15', '2012-01-16', '2012-01-17', '2012-01-18', '2012-01-19', '2012-01-20', '2012-01-21', '2012-01-22', '2012-01-23', '2012-01-24', '2012-01-25', '2012-01-26', '2012-01-27', '2012-01-28', '2012-01-29', '2012-01-30', '2012-01-31', '2012-02-01', '2012-02-02', '2012-02-03', '2012-02-04', '2012-02-05', '2012-02-06', '2012-02-07', '2012-02-08', '2012-02-09', '2012-02-10', '2012-02-11', '2012-02-12', '2012-02-13', '2012-02-14', '2012-02-15', '2012-02-16', '2012-02-17', '2012-02-18', '2012-02-19', '2012-02-20', '2012-02-21', '2012-02-22', '2012-02-23', '2012-02-24', '2012-02-25', '2012-02-26', '2012-02-27', '2012-02-28']
if not os.path.exists('%s_to_%s.pkl'%(dates[0],dates[-1])):
    IM = rm.IonoMap('30d43m17.5ss','21d25m41.9se',dates)
    IM.get_radec_RM(ras,decs)

    datadict = {}
    bm = beamPAPER(150.,nside=16)
    bm_map = np.ma.masked_invalid(hp.orthview(bm,return_projected_map=True))[:,:400]
    plt.close() #negate orthview launching xwindow
    datadict['bm'] = bm_map

    for dayindex in range(IM.RMs.shape[0]):
        datadict[dates[dayindex]] = {}
        #TODO convert keys to JDs? would be neater for everything to be numeric
        for ut in range(IM.RMs.shape[1]):
            map = hp.orthview(IM.RMs[dayindex,ut,:],rot=[0,90],return_projected_map=True)
            map = np.ma.masked_invalid(map)[:,:400] #only get above horizon, mask elsewhere
            datadict[dates[dayindex]][ut] = map
            plt.close() #negate orthview launching xwindow
    output = open('%s_to_%s.pkl'%(dates[0],dates[-1]),'wb')
    pickle.dump(datadict,output)
    output.close()
else: print 'Opening %s_to_%s.pkl'%(dates[0],dates[-1])

pkl = open('%s_to_%s.pkl'%(dates[0],dates[-1]),'rb')
data = pickle.load(pkl)
means,vars,absmeans,absvars = [],[],[],[]
for k in data.keys():
    if k.startswith('20'):
        for i in range(24):
            map = data[k][i]#*data['bm']
            #means.append(np.nanmean(map))
            means.append(map[200,200])
            vars.append(np.nanvar(map))
            absmeans.append(np.nanmean(np.abs(map)))
            absvars.append(np.nanvar(np.abs(map)))
#np.savez('mapdata.npz',means=means,absmeans=absmeans,vars=vars,absvars=absvars)
plt.scatter(vars,means,c=(np.array(range(len(means)))%24)+0.1,alpha=0.7)
plt.colorbar()
plt.xlabel(r'$\sigma^2_{\rm RM(\Omega)}$',size=15)
plt.ylabel(r'$RM(\hat{z})$',size=15)
plt.show()


vvv,mmm = [],[]
for i in range(len(means)):
    hr = (i%24)+2
    if hr<=6 or hr>=18:
        vvv.append(vars[i])
        mmm.append(means[i])
vvv=np.array(vvv)
mmm=np.array(mmm)

counts,xbins,ybins=np.histogram2d(np.sqrt(vvv),mmm,bins=20)
plt.contourf(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()])
plt.colorbar()
plt.xlabel(r'$\sigma_{\rm RM(\Omega)}$',size=15)
plt.ylabel(r'$RM(\hat{z})$',size=15)
plt.show()

import IPython;IPython.embed()

