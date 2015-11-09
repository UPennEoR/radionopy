from __future__ import division, print_function
import os
import sys
import numpy as np
import jdcal
import pylab

def plot_attenuation(rmArray, f, label):
	N = rmArray.shape[0]
	l2 = (0.3 / f) ** 2
	factors = []
	for i, RM in enumerate(rmArray):
		atten = 0
		for j in range(i + 1, N, 1):
			atten += np.cos(2 * (rmArray[j] - RM) * l2)

		factor = (float(N) + (2. * atten)) / pow(float(N), 2.)
		factors.append(factor)

	Pv, v = np.histogram(factors, bins=10, density=True)
	print('Epsilon=', np.mean(factors), '+/-', np.std(factors))
	v = 0.5 * (v[1:] + v[:-1])
	pylab.step(v, Pv, label=label)

	return [np.mean(factors), np.std(factors)]

def plot_attenuation_sim(rmArray, f, daysInSeason=82, Niter=1000, label=''):
	N = daysInSeason
	l2 = (0.3 / f) ** 2
	hist, bins = np.histogram(rm, bins=np.arange(min(rmArray), max(rmArray) + binwidth, binwidth))
	bin_midpoints = bins[:-1] + np.diff(bins) / 2
	cdf = np.cumsum(hist)
	cdf = cdf / cdf[-1]
	values = np.random.rand(Niter)
	value_bins = np.searchsorted(cdf, values)
	random_from_cdf = bin_midpoints[value_bins]

	factors = []

	for i, RM in enumerate(random_from_cdf):
		atten = 0
		for j in range(i + 1, N, 1):
			atten += np.cos(2 * (random_from_cdf[j] - RM) * l2)
		factor = (float(N) + (2. * atten)) / pow(float(N), 2.)
		factors.append(factor)

	Pv, v = np.histogram(factors, bins=10, density=True)
	print('Epsilon (sim) =', np.mean(factors), '+/-', np.std(factors))
	v = 0.5 * (v[1:] + v[:-1])
	pylab.step(v, Pv, label=label)

	return [np.mean(factors), np.std(factors)]

if __name__ == '__main__':
	storage = []
	sast_check = []

	if len(sys.argv[1:]) >= 2:
		for RMfile in sys.argv[1:]:
			hyp_split = os.path.basename(RMfile).split('-')
			year, month, day = int(hyp_split[0]), int(hyp_split[1]), int(hyp_split[2][:2])
			pointing = int(hyp_split[2][2])

			JDtup = jdcal.gcal2jd(year, month, day)
			JD = JDtup[0] + JDtup[1] + 0.5

			#GMST(in hours) = 6.656306 + 0.0657098242 * (JD0-2445700.5) + 1.0027379093 * UT
			#LST = (GMST - (longitude west of Greenwich in degrees) * (24 / 360)) % 24
			#
			#e.g. Santa Cruz is 122.05 degrees west of Greenwich.
			#	  PAPER is -21.428 degrees west of Greenwich

			#print year, month, day, JD, pointing

			UT, TEC, _, RM, eRM = np.loadtxt(RMfile, unpack=True)

			for i, u in enumerate(UT):
				gmst = 6.656306 + 0.0657098242 * (JD - 2445700.5) + 1.0027379093 * u
				lst = (gmst + 21.428 *(24. / 360.)) % 24

				if int(lst) == pointing:
					 a = [JD, RM[i], eRM[i]]
					 storage.append(a)
					 sast_check.append(u + 2)

	LST = np.load('LST1.npz')
	lstdata = LST['data']
	jd1 = lstdata[:, 0]
	rm1 = lstdata[:, 1]
	erm1= lstdata[:, 2]

	LST = np.load('LST4.npz')
	lstdata = LST['data']
	jd4 = lstdata[:, 0]
	rm4 = lstdata[:, 1]
	erm4 = lstdata[:, 2]

	LST = np.load('LST8.npz')
	lstdata = LST['data']
	jd8 = lstdata[:, 0]
	rm8 = lstdata[:, 1]
	erm8 = lstdata[:, 2]

	'''
	print('LST=1: mean RM=', np.mean(rm1), '+/-', np.std(rm1))
	print('LST=4: mean RM=', np.mean(rm4), '+/-', np.std(rm4))
	print('LST=8: mean RM=', np.mean(rm8), '+/-', np.std(rm8))
	'''

	'''
	LST=1 (6pm->8pm SAST): 
	mean RM = 1.5 +/- 0.2
	Epsilon = 0.013 +/- 0.002
	Epsilon (sim, Niter=1000) = 0.012 +/- 0.001

	LST = 4 (7pm->11pm SAST): 
	mean RM = 1.0 +/- 0.2
	Epsilon = 0.014 +/- 0.004
	Epsilon (sim, Niter=1000) = 0.012 +/- 0.001

	LST = 8 (12am->3am SAST): 
	mean RM = 0.7 +/- 0.1
	Epsilon = 0.018 +/- 0.004
	Epsilon (sim, Niter=1000) = 0.013 +/- 0.002
	'''

	binwidth = 0.1

	f, axarr = pylab.subplots(3,sharex=True)

	axarr[0].hist(rm1, bins=np.arange(min(rm1), max(rm1) + binwidth, binwidth))
	axarr[0].text(0.1, 20, 'LST=1h')
	axarr[0].set_xlim(0, 2)
	axarr[0].set_ylim(0, 25)

	axarr[1].hist(rm4, bins=np.arange(min(rm4), max(rm4) + binwidth, binwidth))
	axarr[1].text(0.1, 20, 'LST=4h')
	axarr[1].set_xlim(0, 2)
	axarr[1].set_ylim(0, 25)
	axarr[1].set_ylabel(r'Frequency')

	axarr[2].hist(rm8, bins=np.arange(min(rm8), max(rm8) + binwidth, binwidth))
	axarr[2].text(0.1, 20, 'LST=8h')
	axarr[2].set_xlim(0, 2)
	axarr[2].set_ylim(0 ,25)
	axarr[2].set_xlabel(r'RM (rad m$^{-2}$)')

	pylab.show()
	pylab.close()

	plot_attenuation(rm1, 0.164, label='LST=1')
	plot_attenuation(rm4, 0.164, label='LST=4')
	plot_attenuation(rm8, 0.164, label='LST=8')

	pylab.legend()

	pylab.ylabel('Frequency')
	#pylab.xlabel(r'$\left(N + 2\sum_{i<j}\cos[2(\Phi_i-\Phi_j)\lambda^2]\right)/N^2$')
	pylab.xlabel(r'$\epsilon$',size=15)
	pylab.show()
	pylab.close()

	plot_attenuation_sim(rm1, 0.164,label='LST=1')
	plot_attenuation_sim(rm4, 0.164,label='LST=4')
	plot_attenuation_sim(rm8, 0.164,label='LST=8')

	pylab.ylabel('Frequency')
	#pylab.xlabel(r'$\left(N + 2\sum_{i<j}\cos[2(\Phi_i-\Phi_j)\lambda^2]\right)/N^2$')
	pylab.xlabel(r'$\epsilon_{\rm sim}$',size=15)
	pylab.show()
	pylab.close()
