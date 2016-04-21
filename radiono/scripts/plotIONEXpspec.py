import sys, numpy as np, pylab, healpy as hp

mp = sys.argv[1]

m = hp.read_map(mp)

pylab.hist(m,bins=100)
pylab.xlabel('TEC')
pylab.ylabel('N(TEC)')
pylab.show()

nside = hp.npix2nside(len(m))
lmax = 4*nside
cl = hp.anafast(m,lmax=lmax)
ell = np.arange(len(cl))
p = ell*(ell+1.)*cl/(np.pi*2.)
pylab.loglog(ell,p,'k-')
pylab.xlabel(r'$l$')
pylab.ylabel(r'$C_l$')
pylab.grid()
pylab.show()

import IPython; IPython.embed()
