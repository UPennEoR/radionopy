import numpy as np, healpy as hp, os, sys, optparse, pyfits
from matplotlib import pylab
from datetime import datetime as datetime
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
        #if write2fits: map_cube.write_fits("gsm_cube_%s-%sMHz"%(str(int(flo)), str(int(fhi))))
        M = np.swapaxes(map_cube,0,1) #conform to other map making routines
    return M
    

def mk_fg_cube(onescale=True, pfrac=0.002, flo=100., fhi=200., nbins=203, alo=-2.7, ahi=-2.3, intermediates=True):
    """
    Make a Stokes IQUV cube.
    
    pfrac: fraction Q,U/I
    flo: lowest frequency
    fhi: highest frequenct
    nbins: number of frequency bins
    alo: lowest spectral index
    ahi: highest spectral index
    
    onescale: I'm not feeling smart enough to make an argument interpreter for 
    different "mk_map" cases. If onescale=True, diffuse emission is modelled using a single power law. Otherwise, we use the SCK model.
    
    I'm making the :
        - large assumption that Q and U trace the spectral distribution of 
          diffuse Stokes I power at a fixed polarization fraction.
        - assumption that spectral indicies are randomly distributed on scales > 3deg
        - small assumption that Stokes V power is vanishingly small.
    
    TODO:
        - map Q<->U using polarization angle map
        - maybe move away from assumption that Galactic Sync. is dominant source of polarization?
        - could use a "correlation length" ala Shaw et al. 2014 to get more realistic spectral index distribution
        
    """
    nu = np.linspace(flo,fhi,num=nbins)
    nside=512 #to match GSM nside
    npix = hp.nside2npix(nside) 
    ipix = np.arange(npix)
    
    #Spectral indices -2.7 to -2.3 are the 2sigma range of Rogers & Bowman '08 (in K)
    alpha = np.random.uniform(low=alo,high=ahi,size=npix)
    alpha = hp.smoothing(alpha,fwhm=np.radians(3.))
    
    I = mk_map_GSM(frange=[flo,fhi],nbins=nbins) #use GSM for Stokes I
    #^this call is gonna be a doozy with 203 frequency bins
    
    if onescale:
        Q0 = mk_map_onescale(512) 
        U0 = mk_map_onescale(512) #different realizations of same scaling
        #XXX with a polarization angle map, I could link Q,U instead of having them independent
    else:
        Q0 = mk_map_SCK(512,flo,700,2.4,2.80)
        U0 = mk_map_SCK(512,flo,700,2.4,2.80)
        #XXX as above wrt pol angle, but also this currently assumes we are dominated by Galactic Synchrotron
    
    _Q0,_U0 = Q0 - Q0.min(),U0 - U0.min()
    Q0 = (2*_Q0/_Q0.max()) - 1 #scale to be -1 to 1 
    U0 = (2*_U0/_U0.max()) - 1 #scale to be -1 to 1
    
    if intermediates:
        plot_maps([Q0,U0],titles=['raw Q','raw U'])
        np.savez('rawcube_%sMHz.npz'%str(flo),maps=[I[:,0],Q0,U0])
    
    Qmaps,Umaps,Vmaps = np.zeros((npix,len(nu))),np.zeros((npix,len(nu))),np.zeros((npix,len(nu)))
    
    #If only I could take the log! Then this would be vectorizable
    #stoopid Q and U with their non +ve definition
    for i in ipix:
        Qmaps[i,:] = Q0[i] * np.power(nu/(np.mean([flo,fhi])),alpha[i])
        Umaps[i,:] = U0[i] * np.power(nu/(np.mean([flo,fhi])),alpha[i])
    date = datetime.now().timetuple()
    if intermediates: np.savez('alpha_%i-%i-%i.npz'%(date[0],date[1],date[2]),maps=alpha)
    
    #impose polarization fraction as fraction of sky-average Stokes I power per frequency
    Qmaps *= np.nanmean(I,axis=0)*pfrac
    Umaps *= np.nanmean(I,axis=0)*pfrac
    
    spols = ['I','Q','U','V']
    cube = [I,Qmaps,Umaps,Vmaps]
    for i,m in enumerate(cube):
        N = 'cube_%s_%s-%sMHz.npz'%(spols[i],str(flo),str(fhi))
        print '    Saving %s'%N
        np.savez(N, maps=m)
    
    return np.array(cube)

def propOpp(cube=None,flo=100.,fhi=200.,npznamelist=None):
    """
    Propogate the Q and U components of an IQUV cube through
    the Oppermann et al. 2012 RM map.
    
    The cube must be 4 or 2 by npix by nfreq. If 4, the middle two arrays will be
    assumed to be Q & U IN THAT ORDER. 
    Or if fromnpz=True, then provide an array of npz names (cube length assumptionsremain)
    """
    ## load the maps
    if npznamelist!=None:
        assert(cube==None)
        nNpz = len(npznamelist)
        assert(nNpz == 2 or nNpz == 4)
        
        if nNpz == 2:
            Q = np.load(npznamelist[0])['maps']
            U = np.load(npznamelist[1])['maps']
        else:
            Q = np.load(npznamelist[1])['maps']
            U = np.load(npznamelist[2])['maps']
    elif cube!=None:
        Q = cube[1]
        U = cube[2]
    else:
        raise ImplmentationError('No map information provided.')
    
    ##
    nbins = Q.shape[1]
    nu = np.linspace(flo,fhi,num=nbins)
    lam = 3e8/(nu*1e6)
    lam2 = np.power(lam,2)
    
    d = pyfits.open('opp2012.fits')
    RM = d[3].data.field(0)
    """
    The RM map is nside=128. Everything else is nside=512.
    We're smoothing on scales larger than the pixellization
    this introduces, so no worries. 
    
    RMmap=hp.ud_grade(RM,nside_out=512)
    hp.mollview(RMmap,title='Oppermann map')
    
    RMmap=hp.smoothing(RMmap,fwhm=np.radians(1.))
    
    hp.mollview(RMmap,title='Oppermann map smoothed')
    pylab.show()
    """
    #Downsample RM variance in alm space
    #Upsample it in pixellization
    #Is this kosher?   
    RMmap = hp.alm2map(hp.map2alm(RM,lmax=50),nside=512)
    
    phi = np.outer(RMmap,lam2)
    fara_rot = (Q + 1.j*U)*np.exp(-2.j*phi) #Eq. 9 of Moore et al. 2013
    Qmaps_rot = fara_rot.real
    Umaps_rot = fara_rot.imag
    
    QU = [Qmaps_rot,Umaps_rot]
    
    np.savez('cube_Qrot_%s-%sMHz.npz'%(str(flo),str(fhi)),maps=Qmaps_rot)
    np.savez('cube_Urot_%s-%sMHz.npz'%(str(flo),str(fhi)),maps=Umaps_rot)
    
    return np.array(QU)

def plot_maps(maps,titles=None):
    """
    Plots a list of healpix maps in mollweide projection
    """
    s = int(np.ceil(np.sqrt(float(len(maps)))))
    for i,m in enumerate(maps):
        if titles is not None: hp.mollview(m,title=titles[i],sub=(s,s,i+1))
        else: hp.mollview(m,sub=(s,s,i+1))
    pylab.show()
    
    