import numpy as np
import pylab as plt

nsamp = 100
mu = 1.
sigma = 0.3

def cumAvg(nparray):
    n = np.arange(len(nparray))+1
    return nparray.cumsum()/n

# sanity check
phase = np.random.uniform(low=0.0,high=2.*np.pi,size=nsamp)
randcmplx = np.exp(-1j*phase)

# relative amplitudes
vi = 100.
vq = 1.
vu = 1.


ntries = 100
eps_dist = []
for n in range(ntries):
    rm = np.random.normal(mu,sigma,nsamp)
    lambda2 = 4.
    phi = 2*rm*lambda2
    #phi = np.zeros(nsamp)
    dfm_phs = np.exp(-1j*phi)
    #rotfac = dfm_phs
    rotfac = vi + vq*np.cos(phi) + vu*np.sin(phi)
    eps = 0.
    for i in range(nsamp):
        for j in range(nsamp):
            eps += rotfac[i]*np.conj(rotfac[j])
    eps /= nsamp**2
    eps_dist.append(eps.real)

eps_dist = np.array(eps_dist)
#eps_dist -= vi**2 # subtract off expected I contribution
#eps_dist /= vq**2 # measure relative to original Q power spectrum
print eps_dist.mean()
