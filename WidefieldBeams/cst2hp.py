import numpy as np
import healpy as hp
import pylab as pl

def healpixellize(f_in,theta_in,phi_in,nside,weighted=True):
    """ A brute force method for converting data f sampled at points theta and phi (not on a healpix grid) into a healpix at resolution nside """

    # Input arrays are likely to be rectangular, but this is inconvenient
    f = f_in.flatten()
    theta = theta_in.flatten()
    phi = phi_in.flatten()
    
    pix = hp.ang2pix(nside,theta,phi)

    map = np.zeros(hp.nside2npix(nside))
    hits = np.zeros(hp.nside2npix(nside))
    
    # Default method is gridding based on weighted four-point interpolation to
    # the nearest pixel (from healpy).  Simply finding the nearest pixel is
    # what will happen if weighted is set to False
    if (weighted):
        for i,v in enumerate(f):
            # Find the nearest pixels to the pixel in question
            neighbours,weights = hp.get_neighbours(nside,theta[i],phi[i])
            # Add weighted values to map
            map[neighbours] += v*weights
            # Keep track of weights
            hits[neighbours] += weights
        map = map/hits
        wh_no_hits = np.where(hits == 0)
        map[wh_no_hits[0]] = hp.UNSEEN
    else:    
        for i,v in enumerate(f):
            map[pix[i]] += v
            hits[pix[i]] +=1
        map = map/hits

    return (map,hits)

def cst2hp(cstfile,filetype='rich',column=2,nonorm=False,phase=False):
# It's annoying that the file types are difficul enough to parse that
# I really need completely different sets of parameters

    if filetype == 'eloy':
        skiprows = 3
    if filetype == 'rich':
        skiprows = 2
        
    # Read in the text file
    data = np.loadtxt(cstfile,skiprows=skiprows)
    # Pretty much theta and phi in physics spherical coordinates
    thetad = data[:,0]
    phid = data[:,1]

    if filetype=='eloy':        
        # Eloy's coordinates are really hard to parse in the language of
        # spherical coordinates.  Basically, the first column IS theta, but it
        # wraps to both sides of the pole.  The second column is phi, but you
        # can only figure out what phi is by knowing the signs of both theta
        # and phi.  Here goes.
        tp_pp = np.array(np.where((thetad >= 0) * (phid >=0))[0])
        tp_pm = np.array(np.where((thetad >= 0) * (phid < 0))[0])
        tm_pp = np.array(np.where((thetad < 0) * (phid >= 0))[0])
        tm_pm = np.array(np.where((thetad < 0) * (phid < 0))[0])

        phid_shift = np.zeros(len(phid))
        
        phid_shift[tp_pp] = phid[tp_pp]
        phid_shift[tp_pm] = phid[tp_pm]+360.
        phid_shift[tm_pp] = phid[tm_pp]+180
        phid_shift[tm_pm] = phid[tm_pm]+180

        phid = phid_shift
        thetad = np.abs(thetad)

    # Healpix likes everything in radians
    theta = np.radians(thetad)
    phi = np.radians(phid)

    # For Eloy, third column is E field magnitude, so square and divide by max
    # to get power response
    if filetype=='eloy':
        val = np.power(data[:,column],2)

    # For Rich, third column is dBi    
    if filetype=='rich':
        val = np.power(10.,data[:,column]/10.)

    if not nonorm:
        val = val/val.max()        

    if phase:
        val = data[:,column]
        
    # Original Healpix gridding
    nside=32
    npix = hp.nside2npix(nside)
    
    map, hits = healpixellize(val,theta,phi,nside)

    # Let's decompose this into a_lm and recompose at higher resolution
    if True:
        nside_synth = 128
        lmax = 3*nside-1
        l,m = hp.Alm.getlm(lmax)

        # Do the alm decomposition on dBi map
        alm = hp.map2alm(10.*np.log10(map),lmax=lmax)
        # Reconstitute at high resolution
        bm = hp.alm2map(alm,nside_synth,lmax=lmax,verbose=False)
        bm = np.power(10.,bm/10.)
        bm = bm/bm.max()
        
    return {'data':data,'bm32':map,'bm128':bm}

def theta_cut_at_phi(nside,phi_d_cut,npts=3600):
    """For a HEALpix map of side nside, give the pixel values for a cut that runs -90 < theta < 90 at some phi.  npts should be determinable from nside, but right now I'm lazy"""
    # Ugh.  Usual physics spherical coordinates do not make nice cuts
    # through the beam.  Easiest way to make cuts
    theta = np.linspace(-180,180,npts)
    phi = np.ones_like(theta)
    phi[np.where(theta < 0)] = phi_d_cut+180
    phi[np.where(theta >= 0)] = phi_d_cut

    p = np.radians(phi)
    t = np.radians(np.abs(theta))
    pix = np.array(hp.ang2pix(nside,t,p))
    return {'pix':pix,'theta':theta}
