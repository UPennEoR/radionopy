import fgutils as fgu, healpy as hp
from matplotlib import pylab
#c = fgu.mk_fg_cube(nbins=203,onescale=False)
#qu= fgu.propOpp(npznamelist=['cube_Q_100.0-200.0MHz.npz','cube_U_100.0-200.0MHz.npz'])

flo_fhi_pairs = [[100.,120.,41], [120.5,140.,40], [140.5,160.,41], [160.5,180.,40], [180.5,200.,41]]

for t in flo_fhi_pairs:
    fl,fh,nb = t[0],t[1],t[2]
    c = fgu.mk_fg_cube(flo=fl,fhi=fh,nbins=nb,onescale=False,verbose=True, alpha_map='alpha_2016-3-21.npz',raw_map='rawcube_100.0MHz.npz')
    
    qName = 'cube_%s_%s-%sMHz.npz'%('Q',str(fl),str(fh))
    uName = 'cube_%s_%s-%sMHz.npz'%('U',str(fl),str(fh))
    
    qu= fgu.propOpp(flo=fl,fhi=fh,npznamelist=[qName,uName])

import IPython;IPython.embed()