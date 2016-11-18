from __future__ import division
import numpy as np
import sys
from matplotlib import pyplot as plt

args = sys.argv[1:]
NNN = 1000
cols = ['b','g','r']
freqs = np.linspace(40e6,500e6,num=NNN)

num = 4

stor = np.zeros((num,NNN,len(args)))
rstor= np.zeros((NNN,len(args)))

for c,npz in enumerate(args):
    print npz
    d = np.load(npz)
    rms = d['dhist']
    N = len(rms)
    print np.mean(rms),np.std(rms)
    
    #construct PDF from data
    bw=0.07
    hist,bins=np.histogram(rms,bins=np.arange(min(rms),max(rms)+bw,bw))
    bin_midpoints = bins[:-1] + np.diff(bins)/2
    cdf = np.cumsum(hist)
    cdf = cdf / cdf[-1]
    values = np.random.rand(N)
    value_bins = np.searchsorted(cdf,values)
    #random_from_cdf = bin_midpoints[value_bins]
    rfc = bin_midpoints[value_bins]
    epsarr = [] 
    for n in range(num):
        print n
        grms = np.random.normal(np.mean(rms),scale=np.std(rms),size=(N))
        for k,f in enumerate(freqs):
            lam = 3.e8/f
            lam2 = lam*lam
            rsum,gsum = 0,0
            for i in range(N):
                for j in range(0,i):
                    #gsum+=np.cos( 2*(grms[i]-grms[j])*lam2 )
                    gsum+=np.cos( 2*(rfc[i]-rfc[j])*lam2 ) #drawn from data histograms
                    rsum+=np.cos( 2*(rms[i]-rms[j])*lam2 ) #actual data
            gs,rs = float(N)+(2*gsum),float(N)+(2*rsum)
            geps,reps = gs/float(N*N),rs/float(N*N)
            stor[n,k,c] = geps
            #if f>121e6 and f<131e6: epsarr.append(geps)
            rstor[k,c] = reps
    #print npz,np.mean(epsarr),np.std(epsarr)
    #np.savez(npz.split('.')[0]+'gauss_n%i.npz'%num,res=stor)
    #np.savez(npz.split('.')[0]+'data_n%i.npz'%num,res=stor)
    M,S = np.mean(stor[:,:,c],axis=0),np.std(stor[:,:,c],axis=0)
    plt.errorbar(freqs/1e6,M,yerr=S,label=npz.split('_')[0]+' sim',fmt=cols[c]+'-',ecolor=cols[c])
    plt.plot(freqs/1e6,rstor[:,c],cols[c]+'-',lw=2,alpha=0.7,label=npz.split('_')[0]+' data')
plt.ylabel('Attenuation Factor', size=13)
plt.xlabel('Frequency [MHz]', size=13)
plt.legend()
plt.show()

                
                
            
