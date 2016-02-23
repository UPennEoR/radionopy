# Ugh - he doesn't even have the decency to provide this
import numpy as np
import pylab as plt
import os

#path = '/Users/jaguirre/PyModules/ionFR/'
path = os.path.expanduser('~/radionopy/')

UT, TEC, B, RM, dRM = np.loadtxt(os.path.join(path, 'RM_files/decra.txt'), unpack=True)

plt.figure(1)
plt.clf()
plt.errorbar(UT[0:24], RM[0:24], yerr=dRM[0:24], marker='o', ls='None', color='red')
plt.ylim([0, 3])
plt.xlim([0, 25])

UT, TEC, B, RM, dRM = np.loadtxt(os.path.join(path, 'IonRM_from_ionFR.txt'), unpack=True)
plt.errorbar(UT[0:24], RM[0:24], yerr=dRM[0:24], marker='o', ls='None', color='black')

UT, TEC, B, RM, dRM = np.loadtxt(os.path.join(path, 'IonRM.txt'), unpack=True)
plt.errorbar(UT[0:24], RM[0:24], yerr=dRM[0:24], marker='o', ls='None', color='blue')

#plt.savefig('test.pdf')
plt.show()
