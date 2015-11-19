# Ugh - he doesn't even have the decency to provide this
import numpy as np
import pylab as plt

path = '/Users/jaguirre/PyModules/ionFR/'

UT, TEC, B, RM, dRM = np.loadtxt(path+'IonRM.txt',unpack=True)

plt.clf()
plt.errorbar(UT[0:23],RM[0:23],yerr=dRM[0:23],marker='o',ls='None')
plt.ylim([0,3])
plt.xlim([0,24])
plt.show()
