import numpy as np, matplotlib.pyplot as plt, pickle, sys, os
from radiono import utils as ut
from scipy.stats import norm
import matplotlib.mlab as mlab

if not os.path.exists('lstbintest.npz'):
    print 'Reading pickle...'
    pkl = open('2011-12-06_to_2012-02-28.pkl','rb')
    data = pickle.load(pkl)
    print '    read done.'

    LSTres = 15. #resolution of LST bin in RA degrees. 
    ra_h = np.arange(0.,360.,LSTres) #LST hours
    N = ra_h.shape[0]
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
            try: LSTbins[j,i,:,:] = data[k][pTB[1]%N]
            except IndexError:
                print 'Index error for key %s, ra %i'%(k,int(ra))
                continue

    np.savez('lstbintest.npz',bins=LSTbins)

print 'Opening lstbintest.npz'
LSTbins = np.load('lstbintest.npz')['bins']

print 'Averaging everybody together on the sky'
k=0
stor = np.zeros((LSTbins.shape[0]*LSTbins.shape[1],LSTbins.shape[2],LSTbins.shape[2]))
for i in range(LSTbins.shape[0]):
    for j in range(LSTbins.shape[1]):
        stor[k,:,:] = LSTbins[i,j,:,:]
        k+=1

plt.imshow(np.mean(stor,axis=0))
plt.colorbar()
#plt.show()
plt.close()

plt.imshow(np.std(stor,axis=0))
plt.colorbar()
#plt.show()
plt.close()

f,axarr=plt.subplots(4,4,sharex=True,sharey=True) #XXX neglects low-occupancy LST=16 of PSA32 season


print 'Calculating main lobe RM per LST histograms'

plainHIST = True

for i,ax in enumerate(axarr.ravel()):
    toplot = []
    
    try: a=LSTbins[i,0,0,0]
    except IndexError: continue
    
    for k,entry in enumerate(LSTbins[i,:,150:250,150:250].ravel()):
        if entry!=0.: toplot.append(entry)
    if len(toplot)==0: continue
    bw = 0.05
    B = np.arange(min(toplot),max(toplot)+bw,bw)
    
    (mu,sigma) = norm.fit(toplot)
    print 'LST=%i, mu=%f, sigma=%f'%(i,mu,sigma)
    n, bins, patches = ax.hist(toplot, B)#, normed=1)
    ax.cla()
    
    if plainHIST:
        ax.hist(toplot,B)
        if i%4==0: ax.set_ylabel('Occupancy')
    else:
        y = mlab.normpdf(bins, mu, sigma)#*max(n)
        y = y/y.max()
        y0 = np.zeros_like(y)
        ax.plot(bins,y,'k',lw=2)
        ax.fill_between(bins,y,where=y0<y,facecolor='blue',alpha=0.7)
        if i%4==0: ax.set_ylabel('Normalized Occupancy')
    
    ax.set_title('LST %i'%i)
    if abs(i-len(axarr.ravel())) <= 4: ax.set_xlabel(r'${\rm \phi_{\rm iono}}(\rm \Omega_{ML})  {\rm [rad m^{-2}]}$')
    
#plt.show()
plt.close()


print 'Calculating zj-zi histogram'
f,axarr = plt.subplots(4,4)
bw = 0.01
means,sigs = [],[]
for i,ax in enumerate(axarr.ravel()):
    zen_vals = LSTbins[i,:,200,200]
    zen_vals = zen_vals[zen_vals!=0.]
    if len(zen_vals) < 3: continue
    arr = []
    for k,zi in enumerate(zen_vals):
        for l,zj in enumerate(zen_vals):
            if k==l: continue
            arr.append(zi-zj)
    B = np.arange(min(arr),max(arr)+bw,bw)
    
    (mu,sigma) = norm.fit(arr)
    print 'LST=%i, mu=%f, sigma=%f'%(i,mu,sigma)
    means.append(mu)
    sigs.append(sigma)
    n, bins, patches = ax.hist(arr, B, normed=1)
    y = mlab.normpdf(bins, mu, sigma)
    ax.plot(bins,y,'r--',lw=2)
    
    #ax.hist(arr,B)
    ax.set_title('LST %i'%i)
    ax.set_xlim(-0.3,0.3)
    if i%4==0: ax.set_ylabel('Frequency [%]')
    if abs(i-len(axarr.ravel())) <= 4: ax.set_xlabel(r'${\rm \Delta\phi_{\rm iono}}(\rm \Omega_{zen})  {\rm [rad m^{-2}]}$')
    ax.set_ylim(0,16)
plt.show()
plt.close()

print np.mean(means),np.mean(sigs)
sys.exit()

print 'Plotting time-series in days'

f,axarr = plt.subplots(4,4)
for i,ax in enumerate(axarr.ravel()):
    zen_vals = LSTbins[i,:,200,200]
    zen_vals = zen_vals[zen_vals!=0.]
    ax.plot(range(len(zen_vals)),zen_vals,'b-',lw=2)
    ax.set_title('LST %i'%i)
    ax.set_xlim(0,85)
    ax.set_ylim(0,-0.6)
#plt.show()
plt.close()


