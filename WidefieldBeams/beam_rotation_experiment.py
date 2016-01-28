import numpy as np
import pylab as plt
import healpy as hp
import polarized_beams as pb
reload(pb)
import healpyTools as hpt

# Get the PAPER beams from Rich, which only amplitude info
bradley_paper_file = '/Users/jaguirre/Documents/PAPER/2010_beam/sdipole_05e_eg_ffx_150.txt'
# Get the HERA beams from Rich, which are truly complex
bradley_hera_file = '/Users/jaguirre/PyModules/hera-cst/mdl01/HERA_DISH_paper_feed_cyl36_150mhz'

hera_rich = pb.bradley_hera_polarized_beams(bradley_hera_file)

paper_rich = pb.bradley_paper_polarized_beams(bradley_paper_file)

pbxx = paper_rich['xx']
pbxx = hpt.rotate_healpix_map(pbxx,[0,-120])

#hp.mollview(pbxx)
#hp.mollview(hpt.rotate_healpix_map(pbxx,[0,-120]))
#hp.graticule()
#plt.show()

npix = pbxx.size
nside = hp.npix2nside(npix)
lmax = 3*nside-1

alm_pbxx = hp.map2alm(pbxx,lmax=lmax)
cl_pbxx = hp.alm2cl(alm_pbxx)
l,m= hp.Alm.getlm(lmax)

rot = np.exp(-1j*np.radians(45.)*m)
rotmap = hp.alm2map(alm_pbxx*rot,nside)

#hp.mollview(pbxx)
#hp.mollview(rotmap)
#plt.show()

xydip = pb.xy_ideal_dipole()
# F = f_theta(theta,phi) theta_hat + f_phi phi_hat
# theta_hat = ()xhat + ()yhat + ()zhat
Baltaz = np.array([[np.zeros(npix)],[xydip['xt']],[xydip['xp']]])

hprot = [0,90]

R = hp.Rotator(rot=hprot)
T = R.mat

Bradec = T * Baltaz

Bradec_t = hpt.rotate_healpix_map(Bradec[1,:],hprot)
Bradec_p = hpt.rotate_healpix_map(Bradec[2,:],hprot)

xp,yp,zp = 

