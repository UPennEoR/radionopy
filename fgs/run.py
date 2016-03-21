import fgutils as fgu, healpy as hp, pylab
c = fgu.mk_fg_cube(nbins=203,onescale=False)
qu= fgu.propOpp(npznamelist=['cube_Q_100.0-200.0MHz.npz','cube_U_100.0-200.0MHz.npz'])
import IPython;IPython.embed()