import numpy as np

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
    
    x = np.cos(phi)*np.sin(theta)
    y = np.sin(phi)*np.sin(theta)
    z = np.cos(theta)
    
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