import numpy as np
import sys
from matplotlib import pyplot as plt

cols = ['b','g','r']
for c,npz in enumerate(sys.argv[1:]):
    for f in np.linspace(100e6,200e6,num=203):
        d = np.load(npz)
        rms = d['dhist']
        N = len(rms)
        lam = 3.e8/f
        lam2 = lam*lam
        sum=0
        for i in range(N):
            for j in range(0,i):
                sum+=np.cos( 2*(rms[i]-rms[j])*lam2  )
        s = float(N) + (2*sum)
        eps = s/float(N*N)
        #print npz,f,eps
        plt.plot(f/1e6,eps,cols[c]+'o')
plt.show()
