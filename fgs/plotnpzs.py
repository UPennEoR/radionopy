import numpy as np, healpy as hp, glob
from matplotlib import pylab

freqs = np.linspace(100.,200.,num=205)
Qlist = sorted(glob.glob('cube_Q_*'))
Qrotlist = sorted(glob.glob('cube_Qrot*'))

mQ,mQrot = [],[]

for entry in Qlist:
    print '   reading %s'%entry 
    mQ.append(np.load(entry)['maps'])
for entry in Qrotlist:
    print '   reading %s'%entry
    mQrot.append(np.load(entry)['maps'])

spec = []
for i in range(len(Qlist)): 
    spec.append(mQ[i][1000,:])
spec = np.concatenate(spec)
pylab.plot(freqs,spec)

specrot = []
for i in range(len(Qrotlist)): 
    specrot.append(mQrot[i][1000,:])
specrot = np.concatenate(specrot)   
pylab.plot(freqs,specrot)

pylab.show()

import IPython;IPython.embed()