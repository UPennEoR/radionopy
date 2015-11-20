"""
Simulates the visibilities for a single baseline, given an unpolarized PAPER beam.
Outputs to an npz file.
Authors: James E. Aguirre, Saul A. Kohn
"""
import pylab, random, os, sys, healpy as hp, numpy as np, optparse
from astropy import units as u
from astropy import constants as c
from scipy import integrate
from bm_prms import prms

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


o = optparse.OptionParser()
o.set_usage('vis_sim.py [options] <name of output npz>')
o.set_description(__doc__)
o.add_option('--nside',dest='nside',default=128,help='Healpix resolution.')
o.add_option('-u',dest='u',default=30.,help='magnitude of u vector for baseline (in m).')
o.add_option('-v',dest='v',default=0.,help='magnitude of v vector for baseline (in m).')
o.add_option('-w',dest='w',default=0.,help='magnitude of w vector for baseline (in m).')
o.add_option('--map',dest='map',default=None,help='healpix map to pass overhead. If None, a uniform, 100 K background will be used. If "point", a 500 K point source on a 0 K background will be used.')
o.add_option('--LSTbin',dest='LSTbin',default=False,action='store_true',help='use the PAPER-64 power spectrum LST-binned range and time resolution. Else, do a full transit.')
o.add_option('--not_verbose',dest='noverb',default=False,action='store_true',help='Do *not* print statements')

opts, args = o.parse_args(sys.argv[1:])

#parse some arguments
assert(len(args)>0)
opts.u = float(opts.u)
opts.v = float(opts.v)
opts.w = float(opts.w)

freqs = np.linspace(0.117,0.182,num=131) #aipy likes GHz units. avoiding band edges
nside = opts.nside
npix = hp.nside2npix(nside)

pwd = os.getcwd()
if os.path.exists(pwd+'/XX_beam_maps.npz'):
	beam = np.load('XX_beam_maps.npz')['maps']
	
else:
	#create beams 
	beams = np.zeros((freqs.shape[0],npix))
	###CLEARLY THIS ONLY WORKS FOR LINEAR POLARIZATIONS XX AND YY
	print 'Calculating beams:'
	for i, freq in enumerate(freqs):
		if not opts.noverb: print freq,'GHz'
		bm = prms['beam'](np.array([freq]),nside=nside,lmax=20,mmax=20,deg=7)
		bm.set_params(prms['bm_prms'])
		px = range(hp.nside2npix(nside))
		xyz = hp.pix2vec(nside,px)
		poly = np.array([h.map[px] for h in bm.hmap])
		Axx = np.polyval(poly,freq)
		Axx = np.where(xyz[-1] >= 0, Axx, 0)
		Axx /= Axx.max()
		Axx = Axx*Axx
		beams[i,:] = rotate_hmap(Axx,[21,120]) #[0,0]=north pole, [0,90]=equator, [21,120]=about right for PAPER
	
	np.savez('XX_beam_maps.npz',maps=beams)
	beam = beams


#calculate relevant map parameters
c = 3e8 #m/s
ipix = np.arange(npix)
theta,phi = hp.pix2ang(nside,ipix)

#we care about scales ~21 degrees
lmax=3*nside - 1
l,m = hp.Alm.getlm(lmax)

#frequencies in Hz
nfreq=freqs.shape[0]
nu = np.outer(np.linspace(117e6,182e6,num=nfreq),np.ones(npix))#*u.Hz

#define sky -- completely arbitrary choice of temp
uniform_sky = np.ones(npix)*100.#*u.K

#completely arbitrary choice of noise level XXX UN-USED RIGHT NOW
noise = np.zeros(npix)
for i in range(npix): noise[i] = random.uniform(-100,100)#* u.K

####uniform sky tests and point source tests
if opts.map is None: sky=uniform_sky
elif opts.map is 'point':
	sky = np.zeros(npix)
	#define a point source
	theta0 = np.pi/2.
	phi0 = 0.
	pix0 = hp.ang2pix(nside,theta0,phi0)
	#make it slightly less point-like
	nbs = hp.get_all_neighbours(nside,theta0,phi=phi0)
	sky[pix0]=500#*u.K
	for nb in nbs: sky[nb]=500#*u.K
	sky = rotate_hmap(sky,[180,0])
else: sky = hp.read_map(opts.map)

#promote sky to matrix for frequency axis
sky = np.outer(np.ones(nfreq),sky)*pow(nu/150e6,-0.7)

#decompose sky into alm's
n_alm = len(m)
alm = np.zeros((nfreq,n_alm),dtype='complex128')
print 'Calculating sky a_lm values:'
for i in range(nfreq):
	if not opts.noverb: print nu[i,0]/1e6,'MHz'
	alm[i,:] = hp.map2alm(sky[i,:],lmax=lmax,iter=3)

#calculate fringe factor (true for all freqs)	
s  = np.array(hp.pix2vec(nside,ipix))
b = np.resize(np.repeat(np.array([opts.v,opts.u,opts.w]),npix),[3,npix])#*u.meter
b_dot_s = np.sum(b*s,axis=0)
factor = np.exp(1.j*np.outer(np.ones(nfreq),b_dot_s)*nu/c) #c.c didn't work

#decompose B*fringe into alm's
blm = np.zeros((nfreq,n_alm),dtype='complex128')
print 'Calculating instrument a_lm values:'
for i in range(nfreq):
	if not opts.noverb: print nu[i,0]/1e6,'MHz'
	blm[i,:] = hp.map2alm((beam*factor)[i,:],lmax=lmax,iter=3)

#phasing
if opts.LSTbin: rot_ang = np.linspace(0,5.0592,num=1632)<--LST bin rate
else: rot_ang = np.linspace(-np.pi,np.pi,num=360*4) #<--24 hours

#construct visibilties
n = len(rot_ang)
vis = np.zeros([n,nfreq],dtype='complex128')

print 'Calculating visibilities. Time stamp:'
for i in range(n):
	if not opts.noverb: print i
	rotation = np.outer(np.ones(nfreq),np.exp(-1.j*m*rot_ang[i]))
	vis[i,:] = np.sum(alm*blm*rotation,axis=1)

savefile = sys.argv[1]
print 'Saving visibility to %s...'%savefile
np.savez(savefile,vis=vis)	
