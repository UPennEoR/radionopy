import healpy as hp, pylab, numpy as np, pyfits

#Oppermann+ 2012
#http://wwwmpa.mpa-garching.mpg.de/ift/faraday/2012/index.html
d = pyfits.open('opp2012.fits')

#0th entry is blank
# signal map #the galactic Faraday depth divided by a galactic variance profile
# signal uncertainty 
# faraday depth in rad/m^2
# uncertainty in rad/m^2
# galactic variance profile # the ratio between the Faraday depth and the signal field, in rad/m^2
# reconstructed angular PS of the signal field

contents = ['signal map','signal uncertainty','RM map','dRM map','galactic variance profile','angular PS']

for i in range(0,len(contents)-1):
    hp.mollview(d[i+1].data,title=contents[i])
#pylab.show()
pylab.close()

import IPython; IPython.embed()
