import fgutils as fgu
c = fgu.mk_fg_cube(nbins=10,onescale=False)
qu= fgu.propOpp(npznamelist=['cube_Q_100.0-200.0MHz.npz','cube_U_100.0-200.0MHz.npz'])
import IPython;IPython.embed()