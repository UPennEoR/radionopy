# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 09:08:04 2016

@author: jaguirre
"""

import numpy as np
import pylab as plt
import healpy as hp
import healpyTools as hpt
import polarized_beams as pb
reload(pb)
import rotate_sph_vec as rsv
reload(rsv)

if False:
    # Show that I can rotate zy dipole into xy
    rot = [0,30]
    # Define Healpix coordinates
    nside = 32
    npix = hp.nside2npix(nside)
    ipix = np.arange(npix)
    t,p = hp.pix2ang(nside,ipix)
    # 
    zy = pb.zy_ideal_dipole(nside=nside)
    xy = pb.xy_ideal_dipole(nside=nside)
    
    # Do the rotation
    zy_rxt, zy_rxp = rsv.rotate_vector_field(zy['xt'],zy['xp'],t,p,rot)
    zy_ryt, zy_ryp = rsv.rotate_vector_field(zy['yt'],zy['yp'],t,p,rot)
    
    pb.plot_efield_beams(zy['xt'],zy['xp'],zy['yt'],zy['yp'],figno=1)
    #pb.plot_efield_beams(xy['xt'],xy['xp'],xy['yt'],xy['yp'],figno=2,rot=[0,90])
    pb.plot_efield_beams(zy_rxt,zy_rxp,zy_ryt,zy_ryp,figno=3,rot=[0,0])

# Now rotate the HERA beams to the south
if True:
    cstfile = 'HERA_DISH_paper_feed_cyl36_150mhz'
    if False:
        hera = pb.bradley_hera_polarized_beams(cstfile)
    npix = len(hera['xx'])
    nside = hp.npix2nside(npix)
    ipix = np.arange(npix)
    t,p = hp.pix2ang(nside,ipix)
    
    hera_xta = np.real(hera['xt'])
    hera_xpa = np.real(hera['xp'])
    hera_yta = np.real(hera['yt'])
    hera_ypa = np.real(hera['yp'])
    
    rot = [0,-120]

    pb.plot_efield_beams(hera_xta,hera_xpa,hera_yta,hera_ypa,figno=4,rot=[0,90])    

    hera_rxt, hera_rxp = rsv.rotate_vector_field(hera_xta,hera_xpa,t,p,rot)
    hera_ryt, hera_ryp = rsv.rotate_vector_field(hera_yta,hera_ypa,t,p,rot)
  
    pb.plot_efield_beams(hera_rxt,hera_rxp,hera_ryt,hera_ryp,figno=5,rot=[0,0]) 
  
 #   hera_rxx = np.power(hera_rxt,2) + np.power(hera_rxp,2)


#    hp.mollview(hpt.rotate_healpix_map(hera['xx'],rot=rot),title='XX')
#    hp.mollview(hera_rxta,title='X theta')
#    hp.mollview(hera_rxpa,title='X phi')
#    hp.mollview(hera_rxx,title='sum of squares')

plt.show()

