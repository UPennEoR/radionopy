import numpy as np, healpy as hp, os, sys, optparse
from matplotlib import pylab

"""
Making maps from angular power spectra
"""

def mk_map_onescale(nside,lmin=0,alpha=-2.17,smooth=1.):
    """
    One power law model -> random field realization
    """
    lmax = 3*nside - 1
    ll = np.array(range(lmin,lmax,1))
    Cl = np.power(ll,alpha)
    Cl[0] = 0. #zero mean
    
    Q = hp.synfast(Cl,nside,lmax=lmax)
    Qsm = hp.smoothing(Q,fwhm=np.radians(smooth))
    
    return Qsm

def mk_map_SCK(nside,nu,A,beta,alpha,lmin=0,l_f=1000,nu_f=130,smooth=1.):
    """
    Santos Cooray and Knox (2005) diffuse foreground model -> random field realization
    
    #Fiducial parameters at l_f=1000, nu_f=130 MHz
    #
    # Src (i)                        A (mK^2)  beta  \bar{alpha}  xi
    #
    #Extragalactic point sources     57        1.1   2.07         1.0
    #Extragalactic bremstrahlung     0.014     1.0   2.10         35
    #Galactic synchrotron            700       2.4   2.80         4.0
    #Galactic bremstrahlung          0.088     3.0   2.15         35
    """
    lmax = 3*nside - 1
    ll = np.array(range(lmin,lmax,1))
    Cl = A*np.power(float(l_f)/ll,beta)*np.power(nu_f/nu,2*alpha)
    Cl[0] = 0. #zero mean
    
    Q = hp.synfast(Cl,nside,lmax=lmax)
    Qsm = hp.smoothing(Q,fwhm=np.radians(smooth))
    
    return Qsm

def mk_map_GSM(f=150.,frange=None,nbins=None,write2fits=True):
    """
    Make a Stokes I map/cube using pyGSM
    Can also provide frange=[lo,hi] and nbins to get an image cube (ignores f)
    Both f arguments need to be  in MHz
    """
    from pygsm import GlobalSkyModel
    gsm = GlobalSkyModel()
    
    if frange==None:
        gsm = gsm.generate(f)
        if write2fits: gsm.write_fits("gsm_%sMHz.fits"%str(int(f)))
        M = gsm
    else:
        flo,fhi=frange[0],frange[1]
        try: freqs = np.linspace(flo,fhi,nbins)
        except:
            print 'Did you provide "nbins"?'
            return None
        map_cube = gsm.generate(freqs)
        if write2fits: map_cube.write_fits("gsm_cube_%s-%sMHz"%(str(int(flo)), str(int(fhi))))
        M = map_cube
    return M
    


def plot_maps(maps,titles=None):
    s = int(np.ceil(np.sqrt(float(len(maps)))))
    for i,m in enumerate(maps):
        if titles is not None: hp.mollview(m,title=titles[i],sub=(s,s,i+1))
        else: hp.mollview(m,sub=(s,s,i+1))
    pylab.show()
    
    