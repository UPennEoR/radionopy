import numpy as np
import pylab as plt
import healpy as hp
from astropy import units as u
from astropy import constants as c

import radionopyBeta as ri
reload(ri)
# Not sure if these correspond to the IONEX file below
year = '2012'
month = '01'
day = '26'

# -----------------------
# Step 1: get IONEX file.
# Need to add retrieval of errors and interpolation of position
# -----------------------
IONEXfile = 'CODG0220.12I'
TEC = ri.readIonexTEC(IONEXfile)


if False:
    plt.imshow(TEC['TEC'][0,:,:])
    plt.show()

#AltIon = ri.calcionheight(IONEXfile)
#AltIon = AltIon*1000.0 # km to m
# Read the file once, dipshit
AltIon = TEC['AltIon']

LatObs = np.radians(30.)
AzSource = np.radians(270.)
ZenSource = np.radians(15.)

offLat,offLon,AzPunct,ZenPunct = ri.PuncIonOffset(LatObs,AzSource,ZenSource,AltIon)
