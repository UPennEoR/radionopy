import numpy as np, matplotlib.pyplot as plt, pickle,sys
from radiono import utils as ut

print 'Reading pickle...'
pkl = open('2011-12-06_to_2012-02-28.pkl','rb')
data = pickle.load(pkl)
print '    read done.'

LSTres = 15. #resolution of LST bin in RA degrees. 
ra_h = np.arange(0.,360.,LSTres) #LST hours
ra_rad = np.radians(ra_h)
dec = np.radians(-30.-(43./60.)-(17.5/3600.)) #PAPER/HERA zenith

LSTbins = np.zeros((24,len(data.keys())-1,400,400)) 
#XXX len(data.keys()-1) column will never be totally full. REMEMBER TO PRUNE!
#XXX the 400,400 projection is probably a poor estimate of the RM sky due to horizon projections

for i,k in enumerate(sorted(data.keys())): #indexed by string 'YYYY-MM-DD'
    if not k.startswith('20'): continue #skip past metadata
    
    for j,ra in enumerate(ra_h):
        nT = ut.nextTransit(ut.ion2ephDate(k),ra,dec)
        pTB = ut.parseTransitBasic(nT,SunCheck=True)
        if pTB[2]==True: continue #avoid Sun
        try: LSTbins[j,i,:,:] = data[k][pTB[1]%24]
        except IndexError:
            print 'Index error for key %s, ra %i'%(k,int(ra))
            continue

np.savez('lstbintest.npz',bins=LSTbins)
import IPython;IPython.embed()

LSTbins = np.load('lstbintest.npz')['bins']
f,axarr=plt.subplots(4,5,sharex=True,sharey=True)

for i,ax in enumerate(axarr.ravel()):
    toplot = []
    
    try: a=LSTbins[i,0,0,0]
    except IndexError: continue
    
    for entry in LSTbins[i,:,150:250,150:250].ravel():
        if entry!=0.:
            toplot.append(entry)
    if len(toplot)==0: continue
    bw = 0.05
    B = np.arange(min(toplot),max(toplot)+bw,bw)
    
    ax.hist(toplot,B)
    ax.set_title('LST %i'%i)
plt.show()
    
    
    