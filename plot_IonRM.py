# Ugh - he doesn't even have the decency to provide this
import numpy as np
import pylab as plt
import os

path = '/Users/jaguirre/PyModules/ionFR/'
#path = os.path.expanduser('~/radionopy/')

UT, TEC, B, RM, dRM = np.loadtxt(os.path.join(path, 'IonRM.txt'), unpack=True)

plt.clf()
plt.errorbar(UT[0:24], RM[0:24], yerr=dRM[0:24], marker='o', ls='None')
plt.ylim([0, 3])
plt.xlim([0, 25])
plt.savefig('test.pdf')
plt.show()
