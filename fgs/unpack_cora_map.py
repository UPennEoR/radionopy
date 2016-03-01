import h5py, numpy as np, healpy as hp, pylab

# see: https://github.com/radiocosmology/cora
# created by cora-makesky foreground 128 100.0 200.0 203 cora_map.hdf5

file = 'cora_map.hdf5'
f = h5py.File(file)

I,Q,U,V = f['map'][:,0,:],f['map'][:,1,:],f['map'][:,2,:],f['map'][:,3,:]
spol = ['I','Q','U','V']
for i,pol in enumerate([I,Q,U,V]): np.savez(spol[i]+'.npz',map=pol)