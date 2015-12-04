import sys
# Stupid way to get cst2hp ... need to put this somewhere.
sys.path.append("/Users/jaguirre/PyModules/MyModules/herapy/beam")
import numpy as np
import healpy as hp
import pylab as plt
from cst2hp import cst2hp
from healpyTools import rotate_healpix_map as rotmap
import paper_aipy_beam as bm

# Really not reliable
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def Orth(map):
    hp.orthview(map,half_sky=True)
    return

def Orth90(map):
    hp.orthview(map,half_sky=True,rot=[0,90])
    return

def Show():
    plt.show(())
    return

"""Point the beam at the equator to

1) make use of the fact that we know a z-oriented dipole only has
theta-hat field, so the phi-hat is identically zero.

2) correspondingly, the y-dipole is clearly all phi-hat in the
equatorial plane, and only slowly deviates as you leave the equator

2) Avoid the weird singularities at the pole.  I really want this to
not suffer from numerical defects

"""
apod = False
paper = True
david = False

# This should be wrapped up in a map object
nside = 256
npix = hp.nside2npix(nside)
ipix = np.arange(npix)
theta,phi = hp.pix2ang(nside,ipix)

# Complete WAG
dipole_phase_x = np.exp(-1j*2.*rotmap(phi,[0,90]))
dipole_phase_y = np.exp(-1j*2.*rotmap(rotmap(phi,[0,90]),[90,0]))

D_z = np.sin(theta) # Actual dipole power response for z-oriented dipole
D_y = rotmap(rotmap(D_z,[0,90]),[90,0])
prefix = 'dipole'
if apod:
    prefix = 'apod_dipole'
    apod = rotmap(np.exp(-np.power(theta,2)),[0,-90])
    D_z = D_z * apod
    D_y = D_y * apod
if paper:
    prefix = 'paper'
    #A_paper_xx = (cst2hp('/Users/jaguirre/Documents/PAPER/2010_beam/sdipole_05e_eg_ffx_150.txt',filetype='rich',column=2))['bm128']
    D_z = bm.aipy_hp_beam(0.15,nside,efield=True)
    D_z = rotmap(D_z,[0,-90])
    #D_z = D_z * dipole_phase_x
    D_y = bm.aipy_hp_beam(0.15,nside,efield=True,y=True)
    D_y = rotmap(D_y,[0,-90])
    #D_y = D_y * dipole_phase_y
    
# OK, coordinate system.  "x" here denotes dipole along z, with beam
# directed towards the +x axis.  "y" is actually along the y axis, but
# the beam is still oriented towards +x

#xnorm = np.power(1-np.power(np.sin(theta)*np.cos(phi),2),-0.5)
#Axt = xnorm * np.cos(theta) * np.cos(phi)
#Axp = xnorm * (-1.) * np.sin(phi)
Axt = D_z 
Axp = np.zeros(npix)
ynorm = np.power(1-np.power(np.sin(theta)*np.sin(phi),2),-0.5)
Ayt = D_y * ynorm * np.cos(theta) * np.sin(phi)
if david:
    Ayt = np.zeros(npix)
Ayp = D_y * ynorm * np.cos(phi)

A = np.array([[Axt,Axp],[Ayt,Ayp]])

# Now follow the algebra
AA = np.zeros([4,4,npix],dtype='complex64')
nom_vis = np.zeros([4,4,npix],dtype='complex64')

S = 1./2*np.array([[1.,1,0,0],[0,0,1,1j],[0,0,1,-1j],[1,-1,0,0]])
Sinv = np.array([[1.,0,0,1],[1,0,0,-1],[0,1,1,0],[0,-1j,1j,0]])

for i in range(npix):
    J = A[:,:,i]
    AA[:,:,i] = np.kron(J,np.conj(J))
    nom_vis[:,:,i] = np.dot(Sinv,np.dot(AA[:,:,i],S)) 

#nom_vis = nom_vis.real

