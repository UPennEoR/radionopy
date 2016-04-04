import numpy as np
import math
import healpy as hp

def hs_prod(m1,m2):
    """
    The Hilbert-Schmidt inner product of two matrices m1 and m2:
    <m1,m2> = Tr( m1 m2*)
    """
    return (np.inner(m1,np.conjugate(m2))).trace()

def stokes_project(pol_tensor,x):
    """
    Projects out the desired Stokes parameter from a polarization coherence tensor.
    """
    
    # The polarization tensor basis:
#     sI = np.array([[1,0],[0,1]])/2.
#     sQ = np.array([[1,0],[0,-1]])/2.
#     sU = np.array([[0,1],[1,0]])/2.
#     sV = np.array([[0,-1.j],[1.j,0]])/2.
    
    stokes_parameter = np.empty(pol_tensor.shape[-1],dtype='float64')
    npix = stokes_parameter.shape[0]
    
    if x == 'I':
        sI = np.array([[1,0],[0,1]])/2.
        for i in range(npix):
            stokes_parameter[i] = (hs_prod(sI,pol_tensor[:,:,i])).real
    elif x == 'Q':
        sQ = np.array([[1,0],[0,-1]])/2.
        for i in range(npix):
            stokes_parameter[i] = (hs_prod(sQ,pol_tensor[:,:,i])).real
    elif x == 'U':
        sU = np.array([[0,1],[1,0]])/2.
        for i in range(npix):
            stokes_parameter[i] = (hs_prod(sU,pol_tensor[:,:,i])).real
    elif x == 'V':
        sV = np.array([[0,-1.j],[1.j,0]])/2.
        for i in range(npix):
            stokes_parameter[i] = (hs_prod(sV,pol_tensor[:,:,i])).real
            
    return stokes_parameter 

def rotation_matrix(axis, theta):
    
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    # Taken from the internet: 
    # http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def linear_pol_coord_transform(angle1,angle2,rotation_matrix,Q,U):
    """
    
    Based on the IDL Healpix routine "rotate_coord.pro" by K.M Gorski,
    E. Hivon, and A.J. Banday.
    
    (angle1,angle2) are the angular coordinates of the fields (Q,U) in the
    coordinate system of the basis in which Q and U are defined.
    
    angle1 is the colatitude angle, angle2 is the azimuthal position angle
    
    Written assuming (angle1, angle2) is generated from healpy.pix2ang()
    
    rotation_matrix: the 3x3 rotation matrix that specifies the active rotation
    from the starting coordinate system in which the input (Q,U) are defined
    to the desired coordinate system.
    
    Outputs Qnew,Unew, the linear polarization fields in the basis of the new
    coordinate system.
    
    """
    
    # Will add comments!
    rotmat = np.asarray(rotation_matrix)
    Q = np.asarray(Q)
    U = np.asarray(U)
    
    vec0 = np.array([0,0,1])
    
    x = np.cos(angle2)*np.sin(angle1)
    y = np.sin(angle2)*np.sin(angle1)
    z = np.cos(angle1)
    
    in_uvec = np.stack((x,y,z)) # the cartesian position vectors on the unit sphere.
    
    vec = np.einsum('ab...,b...->a...',rotmat,in_uvec)
    vec0 = np.einsum('ba,b->a',rotmat,vec0)
    
    sin_psi = vec0[1]*in_uvec[0,:] - vec0[0]*in_uvec[1,:]
    cos_psi = -vec0[2] + vec[2,:]*in_uvec[1,:]
    
    norm = sin_psi**2. + cos_psi**2. # this is effectively the determinant of the 
    
    sin_2psi = 2. * sin_psi * cos_psi / norm
    cos_2psi = 1. - 2. * sin_psi**2. / norm
    
    # (Qnew + i Unew) = (Q + i U)e^(2 i psi)
    Qnew = cos_2psi * Q + sin_2psi * U
    Unew = -sin_2psi * Q + cos_2psi * U
    
    return Qnew,Unew

def get_sphr_rotation_matrix(angle1,angle2,rotation_matrix):
    rotmat = np.asarray(rotation_matrix)
    
    vec0 = np.array([0,0,1])
    
    x = np.cos(angle2)*np.sin(angle1)
    y = np.sin(angle2)*np.sin(angle1)
    z = np.cos(angle1)
    
    in_uvec = np.stack((x,y,z)) # the cartesian position vectors on the unit sphere.
    
    vec = np.einsum('ab...,b...->a...',rotmat,in_uvec)
    vec0 = np.einsum('ba,b->a',rotmat,vec0)
    
    sin_psi = vec0[1]*in_uvec[0,:] - vec0[0]*in_uvec[1,:]
    cos_psi = -vec0[2] + vec[2,:]*in_uvec[1,:]
    
    norm = sin_psi**2. + cos_psi**2. # this is effectively the determinant of the 
    
    sin_psi /= np.sqrt(norm)
    cos_psi /= np.sqrt(norm)
    
    return np.array([[cos_psi,sin_psi],[-sin_psi,cos_psi]])

def healpixellize(f_in,theta_in,phi_in,nside,fancy=True):
    """ A dumb method for converting data f sampled at points theta and phi (not on a healpix grid) into a healpix at resolution nside """

    # Input arrays are likely to be rectangular, but this is inconvenient
    f = f_in.flatten()
    theta = theta_in.flatten()
    phi = phi_in.flatten()
    
    pix = hp.ang2pix(nside,theta,phi)

    map = np.zeros(hp.nside2npix(nside))
    hits = np.zeros(hp.nside2npix(nside))
    
    # Simplest gridding is map[pix] = val. This tries to do some
    #averaging Better would be to do some weighting by distance from
    #pixel center or something ...
    if (fancy):
        for i,v in enumerate(f):
            # Find the nearest pixels to the pixel in question
            neighbours,weights = hp.get_interp_weights(nside,theta[i],phi[i])
            # Add weighted values to map
            map[neighbours] += v*weights
            # Keep track of weights
            hits[neighbours] += weights
        map = map/hits
        wh_no_hits = np.where(hits == 0)
        print 'pixels with no hits',wh_no_hits[0].shape
        map[wh_no_hits[0]] = hp.UNSEEN
    else:    
        for i,v in enumerate(f):
            map[pix[i]] += v
            hits[pix[i]] +=1
        map = map/hits

    return map


def hires_healpixellize(map_data,theta,phi,nside_interp,nside_synth):
    
    map = healpixellize(map_data,theta,phi,nside_interp)
    
    lmax = 3*nside_interp-1
#    l,m=hp.Alm.getlm(lmax)
    alm = hp.map2alm(map,lmax=lmax)
    
    return hp.alm2map(alm,nside_synth,lmax=lmax)
     
def test_if_close(array1,array2):
    """
    Returns 0 if the data in the arrays are the same.
    """
    return np.where(np.isclose(array1,array2)==False)[0].shape[0]




