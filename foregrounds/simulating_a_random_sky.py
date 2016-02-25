import numpy as np
import healpy as hp
import pylab as plt
import sys


nside = 64
print 'Pixel resolution w/ nside=%s: %s arcmin'%(str(nside), str(hp.nside2resol(nside,arcmin=True)))
npix = hp.nside2npix(nside)
ipix = np.arange(npix)
lmax = 3*nside-1
l = np.linspace(0,lmax,num=lmax)

"""
#following La Porta & Burigana 2006 interpolation scheme
#assumes C^E_l = C^B_l = C_l = K * l^alpha
cl = []
for il in l:
    if il<6: cl.append(1.)
    if il>=6 and il<=30: cl.append(np.power(il,-2.5))
    if il>30: cl.append(np.power(il,-3.)) 
cl[0] = 0 # Zero mean sky
cl = np.array(cl)
"""

cl = np.power(l,-2.17)
cl[0] = 0.

map = hp.synfast(cl,nside,lmax=lmax)
plt.figure(1)
plt.loglog(l,cl)
hp.mollview(map)
plt.show()


sys.exit()



# Make a pattern of angles with a large quadrupole and some random scatter ...
# Lots of jiggery-pokery here.
l,m = hp.Alm.getlm(lmax)
alm_quad = np.zeros((len(l)),dtype=np.complex128)
wh = np.where((l == 2) * (m == 0))[0]
alm_quad[wh] = 1.+0j
polzn_angle = hp.alm2map(alm_quad,nside)
mn = polzn_angle.min()
mx = polzn_angle.max()
polzn_angle = np.pi/3.*(polzn_angle - mn)/(mx-mn)
polzn_angle += np.random.normal(loc=0.,scale=np.pi/25.,size=npix)
#polzn_angle = polzn_angle%np.pi

Q = hp.synfast(cl,nside,lmax=lmax)
#U = Q * np.tan(polzn_angle) # Might be a factor of 2 here ...
U = hp.synfast(cl,nside,lmax=lmax)

# Could generate this with "correlation length" using synfast
alpha = np.random.uniform(low=-1.0,high=-0.7,size=npix)
alpha = hp.smoothing(alpha,fwhm=np.radians(3.))

# Frequencies and wavelengths
nu = np.linspace(100,200,num=203)
lam = 3e8/(nu*1e6)
lam2 = np.power(lam,2)

Qmaps = np.zeros((npix,len(nu)))
Umaps = np.zeros((npix,len(nu)))

for i in ipix:
    Qmaps[i,:] = Q[i] * np.power(nu/150.,alpha[i])
    Umaps[i,:] = U[i] * np.power(nu/150.,alpha[i])

# Need an RM map
RMmap = np.ones_like(Q)*10.
phi = np.outer(RMmap,lam2)

# Faraday rotate
fara_rot = (Qmaps + 1j * Umaps)*np.exp(-1j*phi)
Qmaps_rot = fara_rot.real
Umaps_rot = fara_rot.imag

    
if False:
    plt.figure(1)
    plt.loglog(l,cl)
    hp.mollview(map)
    plt.show()

plt.plot(nu,Qmaps_rot[1000,:])
plt.plot(nu,Umaps_rot[1000,:])
plt.show()
