import pylab, numpy as np

#What is the attenuation of a single V_Q (t,nu) LST bin
#due to ionospheric RM fluctuations?

Ndays = 80
freq = 0.15 #GHz
lamd = 3.e8/(0.15e9)

Nobs = np.linspace(0,Ndays,Ndays)

#I don't know what the center of the RM distribution is
#eyeballing the LST=4 distribution for Moore+'16 RM~0.9
#sigmaRM~0.2
RMdist = np.random.normal(0.9,0.2,Ndays)

noise_curve = pow(Nobs,-0.5)

fiducial_Qvis = 1.+0.j

atten_curve = np.zeros(Ndays)
for i,rm in enumerate(RMdist):
    factor = np.exp(-1.j*rm/pow(lamd,2.))
    if i == 0: atten_curve[i] = fiducial_Qvis*factor
    else:
        temp = np.append(atten_curve[0:i],fiducial_Qvis*factor)
        atten_curve[i] = np.mean(temp)

pylab.plot(Ndays,noise_curve,'k-',label='Pure noise')
pylab.plot(Ndays,np.abs(atten_curve),'b--',label='atten')
pylab.legend()
pylab.show()


        
