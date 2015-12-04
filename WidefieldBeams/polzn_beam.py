import numpy as np
import healpy as hp
import pylab as plt
from cst2hp import cst2hp
from rotate_healpix_map import rotate_healpix_map as rotmap

# Super brute force.  Take a vector function on the sphere
# (theta_hat,phi_hat maps), and calculate the components in the new
# coordinates system defined by rot.
def rotate_vec_func(tmap,pmap,rot):
    npix = len(tmap)
    nside = hp.npix2nside(npix)
    ipix = np.arange(npix)
    theta,phi = hp.pix2ang(nside,ipix)
    # Canonical definitions
    theta_hat = np.array([[np.cos(theta)*np.cos(phi)],[np.cos(theta)*np.sin(phi)],[-np.sin(theta)]]).reshape(3,npix)
    phi_hat = np.array([[-np.sin(phi)],[np.cos(phi)],[np.zeros(npix)]]).reshape(3,npix)
    r = hp.Rotator(rot=rot)
    theta_rot = r(theta_hat)
    phi_rot = r(phi_hat)

    rotmat = np.zeros([2,2,npix])
    rotmat[0,0,:] = (theta_rot*theta_hat).sum(axis=0)
    rotmat[0,1,:] = (theta_rot*phi_hat).sum(axis=0)
    rotmat[1,0,:] = (phi_rot*theta_hat).sum(axis=0)
    rotmat[1,1,:] = (phi_rot*phi_hat).sum(axis=0)

    # Create the new field
    tmap_new = rotmat[0,0]*tmap + rotmat[0,1]*pmap
    pmap_new = rotmat[1,0]*tmap + rotmat[1,1]*pmap

    return {'tmap':tmap_new,'pmap':pmap_new,'rotmat':rotmat}

def Orth(map):
    hp.orthview(map,rot=[0,90],half_sky=True)
    return

nside = 32 #256
npix = hp.nside2npix(nside)
ipix = np.arange(npix)

theta,phi = hp.pix2ang(nside,ipix)

# Get the PAPER beams.  Somethin' weird about these
#A_paper_xx = cst2hp('/Users/jaguirre/Documents/PAPER/2010_beam/sdipole_05e_eg_ffx_150.txt',filetype='rich',column=2)
#A_paper_xt = cst2hp('/Users/jaguirre/Documents/PAPER/2010_beam/sdipole_05e_eg_ffx_150.txt',filetype='rich',column=3)
#A_paper_xp = cst2hp('/Users/jaguirre/Documents/PAPER/2010_beam/sdipole_05e_eg_ffx_150.txt',filetype='rich',column=4)

# A half-wave dipole along z has a pure theta-hat E-field radiation pattern
D_z = np.exp(-np.power(theta,2))
#np.cos(np.pi/2.*np.cos(theta))/np.sin(theta)
# What does this thing look like if I put it along x?
# Rotate from the coordinates of the original problem to the new one
#D_x = D_z #rotmap(D_z,[0,90])
#exiD_y = D_z #rotmap(D_x,[90,0])
# Now, re-express theta-hat in the old system as a linear combination
# of theta and phi in the new one ...
PAPER=False

if PAPER:
    Axt = np.sqrt(A_paper_xt['bm32']/A_paper_xx['bm32'].max()) 
    Axp = np.sqrt(A_paper_xp['bm32']/A_paper_xx['bm32'].max()) 
    Ayt = rotmap(Axt,[90,0])
    Ayp = rotmap(Axp,[90,0])
else:
    xnorm = np.power(1-np.power(np.sin(theta)*np.cos(phi),2),-0.5)
    ynorm = np.power(1-np.power(np.sin(theta)*np.sin(phi),2),-0.5)
    #Axt = D_x * xnorm * np.cos(theta)*np.cos(phi)
    #Axp = D_x * xnorm * -1.*np.sin(phi)
    #Ayt = D_y * ynorm * np.cos(theta)*np.sin(phi)
    #Ayp = D_y * ynorm * np.cos(phi)

    Axt = np.abs(D_z * np.cos(phi))
    Axp = np.abs(D_z * np.sin(phi))
    Ayt = np.abs(D_z * (-np.sin(phi)))
    Ayp = np.abs(D_z * np.cos(phi))

    
Ax_dot_Ay = Axt*Ayt + Axp*Ayp

A = np.array([[Axt,Axp],[Ayt,Ayp]])

# Let's just construct the 4x4 matrix from the algebra, and not use
# np.kron or anything else

A_UU = Axp*Ayt + Axt*Ayp
A_QQ = Axt**2 - Axp**2 - Ayt**2 + Ayp**2

#Orth(A_UU)

#plt.show()

# Now, I basically want to form the 2x2 matrix
# [[Axt,Axp],[Ayt,Ayp]] and outer product it with itself to get the 16 beams
# Brute force ...

AA = np.zeros([4,4,npix],dtype='complex64')
nom_vis = np.zeros([4,4,npix],dtype='complex64')

S = 1./2*np.array([[1.,1,0,0],[0,0,1,1j],[0,0,1,-1j],[1,-1,0,0]])
Sinv = np.array([[1.,0,0,1],[1,0,0,-1],[0,1,1,0],[0,-1j,1j,0]])

for i in range(npix):
    J = A[:,:,i]
    AA[:,:,i] = np.kron(J,J)
    nom_vis[:,:,i] = np.dot(Sinv,np.dot(AA[:,:,i],S)) 

nom_vis = nom_vis.real
    
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
    
beam_labels = np.array([[r'A_x^{\theta}',r'A_x^\phi'],[r'A_y^{\theta}',r'A_y^\phi']])
plt.figure(1)
for a in range(2):
    for b in range(2):
        hp.orthview(A[a,b,:],rot=[0,90],half_sky=True,sub=(2,2,2*a+b+1),min=-1,max=1)
        plt.title(beam_labels[a,b])
        hp.graticule()

plt.figure(2)
for a in range(4):
    for b in range(4):
        hp.orthview(nom_vis[a,b,:],rot=[0,90],half_sky=True,title='',sub=(4,4,4*a+b+1),min=0,max=1)
        hp.graticule()

plt.show()
