import numpy as np
import sys
from matplotlib import pyplot as plt

cols = ['b','g','r']

l4,l6,l8,g4,g6,g8=[],[],[],[],[],[]
llist = [l4,l6,l8]
glist = [g4,g6,g8]

freqs = np.linspace(100e6,200e6,num=203)

for c,npz in enumerate(sys.argv[1:]):
    print npz
    d = np.load(npz)
    rms = d['dhist']
    N = len(rms)
    print np.mean(rms),np.std(rms)
    grms = np.random.normal(np.mean(rms),scale=np.std(rms),size=(N))
    
    for f in freqs:
        lam = 3.e8/f
        lam2 = lam*lam
        gsum,sum=0,0
            
        for i in range(N):
            for j in range(0,i):
                sum+=np.cos( 2*(rms[i]-rms[j])*lam2 )
                gsum+=np.cos( 2*(grms[i]-grms[j])*lam2 )
        s = float(N) + (2*sum)
        gs= float(N) + (2*gsum)
        eps = s/float(N*N)
        geps= gs/float(N*N)
        llist[c].append(eps)
        glist[c].append(geps)

for n in range(3):
    plt.plot(freqs/1e6,llist[n],cols[n]+'-')
    plt.plot(freqs/1e6,glist[n],cols[n]+'--')
plt.show()
