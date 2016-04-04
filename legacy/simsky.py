import numpy as np, pylab, sys, healpy as hp, pyfits
from scipy.interpolate import interp1d

nside = 128 #int(sys.argv[1])
print 'Pixel resolution w/ nside=%s: %s arcmin'%(str(nside), str(hp.nside2resol(nside,arcmin=True)))
npix = hp.nside2npix(nside)
ipix = np.arange(npix)
lmax = 3*nside-1

#Bernardi unpublished APSs (different RMs, take a median)
l_B = [150,200,250,300,350,400,450,500,550,700,800]
l2Cl_2pi_1 = [2e3,3e3,3.5e3,4e3,4.5e3,5e3,5.5e3,6e3,5.5e3,4e3,3.5e3]
l2Cl_2pi_2 = [4e3,6.5e3,8e3,9e3,1e4,1.15e4,1.15e4,1e4,1e4,8e3,5e3]
l2Cl_2pi_3 = [2.5e4,3.5e4,4e4,4e4,4e4,3.5e4,3.5e4,3.5e4,3e4,2e4,1.5e4]
l2Cl_2pi_4 = [7.5e4,8e4,7e4,6e4,6e4,6e4,5.5e4,4.5e4,4e4,2.5e4,2e4]

l2Cl_2pi_med = np.median(np.array([l2Cl_2pi_1,l2Cl_2pi_2,l2Cl_2pi_3,l2Cl_2pi_4]),axis=0)
Cl_med = l2Cl_2pi_med*np.pi*2./np.power(l_B,2.)

l_lo=range(0,11,1)
l_hi=range(900,2800,100)
l_B = l_lo+l_B+l_hi
"""
Cl_lo = []
Cl_hi = []
#La Porta & Burigana Table 2, extrap. to 150 MHz. Valid??
a1 = -0.003*151 + 0.63
a2 = -0.001*151 - 1.29 #<-- does not give sensible results
for L in l_lo: 
    if L<2: Cl_lo.append(0.) #???
    elif L>=2 and L<=11: Cl_lo.append(11-L/10.)#Cl_lo.append(np.power(L,a1))
    #else: Cl_lo.append(3*(11/L))#Cl_lo.append(np.power(10,3.5)*np.power(L,a2))
    
for L in l_hi: Cl_hi.append(0.1*np.power(L/700.,-1.65)) #Bernardi et al. 2009 (rescaled for now)
Cl = np.concatenate((Cl_lo,Cl_med,Cl_hi))
"""
Cl = np.power(np.array(l_B),-1.65)
#f2 = interp1d(l_B,Cl,kind='cubic')
f2 = interp1d(l_B,Cl,kind='linear')
l_new = np.linspace(l_B[0],l_B[-1],num=lmax,endpoint=True)
Cl_new = f2(l_new)
pylab.loglog(l_B,Cl,'o')
pylab.loglog(l_new,Cl_new)
pylab.show()

#sys.exit()

### generate Q and U from C_l

Q,U = hp.synfast(Cl_new,nside,lmax=lmax),hp.synfast(Cl_new,nside,lmax=lmax)
Q,U = hp.smoothing(Q,fwhm=np.radians(1.)),hp.smoothing(U,fwhm=np.radians(1.))

### generate spectral indicies
# Could generate this with "correlation length" using synfast
alpha = np.random.uniform(low=-1.0,high=-0.7,size=npix)
alpha = hp.smoothing(alpha,fwhm=np.radians(3.))

# Frequencies and wavelengths
nu = np.linspace(100,200,num=203)
lam = 3e8/(nu*1e6)
lam2 = np.power(lam,2)

Qmaps = np.zeros((npix,len(nu)))
Umaps = np.zeros((npix,len(nu)))

#If only I could take the log! Then this would be linearizable
#stoopid Q and U with their non +ve definition
for i in ipix:
    Qmaps[i,:] = Q[i] * np.power(nu/150.,alpha[i])
    Umaps[i,:] = U[i] * np.power(nu/150.,alpha[i])

# Need an RM map -- using Opperman et al. (2012)
# http://wwwmpa.mpa-garching.mpg.de/ift/faraday/2012/index.html
d = pyfits.open('opp2012.fits')
RMmap = d[3].data.field(0)
phi = np.outer(RMmap,lam2)

# Faraday rotate
fara_rot = (Qmaps + 1j * Umaps)*np.exp(-1j*phi)
Qmaps_rot = fara_rot.real
Umaps_rot = fara_rot.imag


hp.mollview(Q,title='Q')
hp.mollview(U,title='U')
hp.mollview(alpha,title='alpha')
hp.mollview(RMmap,title='Opperman 2012')


hp.mollview(Qmaps_rot[:,100],title='Q rot 150 MHz')
hp.mollview(Umaps_rot[:,100],title='U rot 150 MHz')
pylab.show()


import IPython; IPython.embed()