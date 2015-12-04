path = '/Users/jaguirre/PyModules/radpolpy/'

# Plot up the individual beams
beam_labels = np.array([[r'A_x^{\theta}',r'A_x^\phi'],[r'A_y^{\theta}',r'A_y^\phi']])

minv = np.array([[0,0],[-0.3,0]])
maxv = np.array([[1,1],[0.3,1]])
        
plt.figure(1)
for a in range(2):
    for b in range(2):
        hp.orthview(np.abs(A[a,b,:]),half_sky=True,sub=(2,2,2*a+b+1),\
                    min=minv[a,b],max=maxv[a,b])
        #plt.title(beam_labels[a,b])
        hp.graticule()

plt.savefig(path+prefix+'_beams.pdf')

# Plot up the "polarization beams"
plt.figure(2)
for a in range(4):
    for b in range(4):
        hp.orthview(np.log10(np.abs(nom_vis[a,b,:])),half_sky=True,title='',sub=(4,4,4*a+b+1),min=-2,max=0)
        hp.graticule()
plt.savefig(path+prefix+'_leakage_beams.pdf')
        
#plt.figure(3)
#moore = Axt**2 + Axp**2 - (Ayt**2 + Ayp**2)
#Orth(moore)
#hp.graticule()

plt.show(())

