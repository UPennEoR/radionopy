import numpy as np
import healpy as hp
import healpyTools as hpt
import pylab as plt

def visibility_leakage_beams(Axt,Axp,Ayt,Ayp):
    npix = Axt.size
    # Define antenna Jones matrix
    A = np.array([[Axt,Axp],[Ayt,Ayp]])

    # Now follow the algebra
    AA = np.zeros([4,4,npix],dtype='complex64')
    stokes_vis = np.zeros([4,4,npix],dtype='complex64')

    S = 1./2*np.array([[1.,1,0,0],[0,0,1,1j],[0,0,1,-1j],[1,-1,0,0]])
    Sinv = np.array([[1.,0,0,1],[1,0,0,-1],[0,1,1,0],[0,-1j,1j,0]])

    for i in range(npix):
        J = A[:,:,i]
        AA[:,:,i] = np.kron(J,np.conj(J))
        stokes_vis[:,:,i] = np.dot(Sinv,np.dot(AA[:,:,i],S))

    return stokes_vis

def dbi2pknorm(dbi):
    pknorm = np.power(10.,dbi/10.)
    pknorm = pknorm/pknorm.max()
    return pknorm

# Mother fucking fucker.  I'm basically going to have to re-write
# cst2hp to properly deal with the polarized beams.  Just like
# everyting else with polarization.
def bradley_paper_polarized_beams(cstfile):
    # This only works on files like /Users/jaguirre/Documents/PAPER/2010_beam/sdipole_05e_eg_ffx_150.txt
    data = np.loadtxt(cstfile,skiprows=2)
    theta = np.radians(data[:,0])
    phi = np.radians(data[:,1])
    # Pixellize in linear space
    dbi = dbi2pknorm(data[:,2])
    dbi_t = dbi2pknorm(data[:,3])
    dbi_p = dbi2pknorm(data[:,4])

    nside=32
    npix = hp.nside2npix(nside)

    xx = hpt.healpixellize(dbi,theta,phi,nside)
    # These are the power beams for each basis vector
    xt = np.sqrt(hpt.healpixellize(dbi_t,theta,phi,nside))
    xp = np.sqrt(hpt.healpixellize(dbi_p,theta,phi,nside))

    yy = hpt.rotate_healpix_map(xx.copy(),[90,0]) # This is a simple rotation, surely
    yt = hpt.rotate_healpix_map(xt.copy(),[90,0])
    yp = hpt.rotate_healpix_map(xp.copy(),[90,0])
    
    return {'xx':xx,'yy':yy,'xt':xt,'xp':xp,'yt':yt,'yp':yp}#'data':data}

def bradley_hera_polarized_beams(cstfile):
    nside=32
    npix = hp.nside2npix(nside)

    data = np.loadtxt(cstfile+'_X.txt',skiprows=2)
    theta = np.radians(data[:,0])
    phi = np.radians(data[:,1])
    # Pixellize in linear space
    dbi = dbi2pknorm(data[:,2])
    dbi_t = dbi2pknorm(data[:,3])
    phase_t = np.radians(data[:,4])
    dbi_p = dbi2pknorm(data[:,5])
    phase_p = np.radians(data[:,6])
    
    xx = hpt.healpixellize(dbi,theta,phi,nside).real
    # These are the power beams for each basis vector
    xt = np.sqrt(hpt.healpixellize(dbi_t,theta,phi,nside))
    xt_phase = hpt.healpixellize(phase_t,theta,phi,nside)
    xp = np.sqrt(hpt.healpixellize(dbi_p,theta,phi,nside))
    xp_phase = hpt.healpixellize(phase_p,theta,phi,nside)
    xt = xt * np.exp(-1j*xt_phase)
    xp = xp * np.exp(-1j*xp_phase)

    data = np.loadtxt(cstfile+'_Y.txt',skiprows=2)
    theta = np.radians(data[:,0])
    phi = np.radians(data[:,1])
    # Pixellize in linear space
    dbi = dbi2pknorm(data[:,2])
    dbi_t = dbi2pknorm(data[:,3])
    phase_t = np.radians(data[:,4])
    dbi_p = dbi2pknorm(data[:,5])
    phase_p = np.radians(data[:,6])

    yy = hpt.healpixellize(dbi,theta,phi,nside).real
    # These are the power beams for each basis vector
    yt = np.sqrt(hpt.healpixellize(dbi_t,theta,phi,nside))
    yt_phase = hpt.healpixellize(phase_t,theta,phi,nside)
    yp = np.sqrt(hpt.healpixellize(dbi_p,theta,phi,nside))
    yp_phase = hpt.healpixellize(phase_p,theta,phi,nside)
    yt = yt * np.exp(-1j*yt_phase)
    yp = yp * np.exp(-1j*yp_phase)
    
    return {'xx':xx,'yy':yy,'xt':xt,'xp':xp,'yt':yt,'yp':yp}

def xy_ideal_dipole():
    nside=32
    npix = hp.nside2npix(nside)
    ipix = np.arange(npix)
    theta,phi = hp.pix2ang(nside,ipix)
    st = np.sin(theta)
    ct = np.cos(theta)
    sp = np.sin(phi)
    cp = np.cos(phi)

    xnorm = 1./np.power((1-np.power(st*cp,2)),0.5)
    ynorm = 1./np.power((1-np.power(st*sp,2)),0.5)
    
    pxt = xnorm * ct * cp
    pxp = xnorm * (-sp)
    pyt = ynorm * ct * sp 
    pyp = ynorm * cp

    return {'xt':pxt,'xp':pxp,'yt':pyt,'yp':pyp}

def zy_ideal_dipole():
    nside=32
    npix = hp.nside2npix(nside)
    ipix = np.arange(npix)
    theta,phi = hp.pix2ang(nside,ipix)
    st = np.sin(theta)
    ct = np.cos(theta)
    sp = np.sin(phi)
    cp = np.cos(phi)

    xnorm = 1./np.power((1-np.power(st*cp,2)),0.5)
    ynorm = 1./np.power((1-np.power(st*sp,2)),0.5)

    # Maybe?
    pxt = np.ones(npix) #xnorm * ct * cp
    pxp = np.zeros(npix) #xnorm * (-sp)
    pyt = ynorm * ct * sp 
    pyp = ynorm * cp

    return {'xt':pxt,'xp':pxp,'yt':pyt,'yp':pyp}

def plot_efield_beams(xt,xp,yt,yp,rot=[0,0],minv = [-1,-1,-1,-1],maxv=[1,1,1,1],figno=1):
    
    A = [xt,xp,yt,yp]
    plt.figure(1)
    plt.clf()
    for a in range(4):
        hp.orthview(A[a],half_sky=True,sub=(2,2,a+1),min=minv[a],max=maxv[a],rot=rot)
        #plt.title(beam_labels[a,b])
        hp.graticule()
#    plt.show()
    return
        
#plt.savefig(path+prefix+'_beams.pdf')

def plot_leakage_beams(nom_vis,rot=[0,0],mn=-2,mx=0,log=True,figno=2):
    if log:
        nom_vis_to_plot = np.log10(np.abs(nom_vis.copy()))
    else:
        nom_vis_to_plot = nom_vis.copy()
    fig = plt.figure(figno)
    plt.clf()
    for a in range(4):
        for b in range(4):
            if (a==b) and (not log):
                mn_plot = 0
                mx_plot = 1
            else:
                mn_plot = mn
                mx_plot = mx
            print mn_plot, mx_plot
            hp.orthview(nom_vis_to_plot[a,b,:],half_sky=True,title='',sub=(4,4,4*a+b+1),min=mn_plot,max=mx_plot,rot=rot)#,margins=[0.1,0.9,0.9,0.1])
        hp.graticule()
    #fig.tight_layout()
#    plt.show()

    return
