import pylab, aipy,sys
import numpy as np

f = np.load(sys.argv[1])
vis = f['vis']

T,N = vis.shape

W = np.ones(N)
#freqs = np.linspace(0.1,0.2,num=N)
freqs = np.linspace(0.112,0.187,num=N)
delays = np.fft.fftfreq(N,d=np.diff(freqs)[0])
weights = aipy.dsp.gen_window(N,'blackman-harris')

delay_vis = np.zeros_like(vis)
fring_vis = np.zeros_like(vis)

for t in range(T):
	_d,info = aipy.deconv.clean(np.fft.ifft(vis[t,:]*weights),np.fft.ifft(W),tol=1e-3)
	dspec = np.fft.fftshift(np.abs(_d+info['res']))
	#dspec = np.fft.fftshift(_d+info['res'])
	delay_vis[t,:] = dspec

W = np.ones(T)
times = range(T)
rates = np.fft.fftfreq(T,d=np.diff(freqs)[0])
weights = aipy.dsp.gen_window(T,'blackman-harris')


for f in range(N):
	_d,info = aipy.deconv.clean(np.fft.ifft(vis[:,f]*weights),np.fft.ifft(W),tol=1e-3)
	fspec = np.fft.fftshift(np.abs(_d+info['res']))
	fring_vis[:,f] = fspec

f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = pylab.subplots(3, 2, sharex='col')#, sharex='col', sharey='row')

ax1.imshow(vis.real,aspect='auto',interpolation='None')
ax1.set_title(r'Re(V$_{ij}$)')

ax2.imshow(vis.imag,aspect='auto',interpolation='None')
ax2.set_title(r'Im(V$_{ij}$)')

ax3.imshow(np.angle(vis),aspect='auto',interpolation='None')
ax3.set_title(r'arg(V$_{ij}$)')

ax4.imshow(np.absolute(vis),aspect='auto',interpolation='None')
ax4.set_title(r'|V$_{ij}$|')

ax5.imshow(np.absolute(delay_vis),aspect='auto',interpolation='None')
ax5.set_title(r'D(V$_{ij}$)')

ax6.imshow(fring_vis.real,aspect='auto',interpolation='None')
ax6.set_title(r'F(V$_{ij}$)')

pylab.show()

ddr_vis = np.zeros_like(vis)

for f in range(N):
	_d,info = aipy.deconv.clean(np.fft.ifft(delay_vis[:,f]*weights),np.fft.ifft(W),tol=1e-3)
	ddrspec = np.fft.fftshift(np.abs(_d+info['res']))
	ddr_vis[:,f] = ddrspec

pylab.imshow(ddr_vis.real,aspect='auto',interpolation='None')
pylab.colorbar()
pylab.show()
pylab.close()