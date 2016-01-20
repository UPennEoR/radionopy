import numpy as np
import pylab as plt
import polarized_beams as pb
reload(pb)
import healpyTools as hpt

# Get the PAPER beams from Rich, which only amplitude info
bradley_paper_file = '/Users/jaguirre/Documents/PAPER/2010_beam/sdipole_05e_eg_ffx_150.txt'
# Get the HERA beams from Rich, which are truly complex
bradley_hera_file = '/Users/jaguirre/PyModules/hera-cst/mdl01/HERA_DISH_paper_feed_cyl36_150mhz'

hera_rich = pb.bradley_hera_polarized_beams(bradley_hera_file)

paper_rich = pb.bradley_paper_polarized_beams(bradley_paper_file)
l_paper_rich = pb.visibility_leakage_beams(paper_rich['xt'],paper_rich['xp'],paper_rich['yt'],paper_rich['yp'])
pb.plot_leakage_beams(l_paper_rich,rot=[0,90],mn=-2,mx=0,log=True,figno=1)

# The sqrt of xx times the ideal dipole response isn't a bad
# approximation for either theta or phi (few percent difference)
# But crucially, the ideal dipole is *not* positive definite
pxy = pb.xy_ideal_dipole()
pdxt = np.sqrt(paper_rich['xx'])*pxy['xt']
pdxp = np.sqrt(paper_rich['xx'])*pxy['xp']
pdyt = np.sqrt(paper_rich['yy'])*pxy['yt']
pdyp = np.sqrt(paper_rich['yy'])*pxy['yp']
l_paper_ideal = pb.visibility_leakage_beams(pdxt,pdxp,pdyt,pdyp)
pb.plot_leakage_beams(l_paper_ideal,rot=[0,90],mn=-2,mx=0,log=True,figno=2)

# Let's try same thing, but at the equator
pzy = pb.zy_ideal_dipole()
pdxt_zy = hpt.rotate_healpix_map(np.sqrt(paper_rich['xx']),[0,-90])*pzy['xt']
pdxp_zy = hpt.rotate_healpix_map(np.sqrt(paper_rich['xx']),[0,-90])*pzy['xp']
pdyt_zy = hpt.rotate_healpix_map(np.sqrt(paper_rich['yy']),[0,-90])*pzy['yt']
pdyp_zy = hpt.rotate_healpix_map(np.sqrt(paper_rich['yy']),[0,-90])*pzy['yp']
l_paper_ideal_zy = pb.visibility_leakage_beams(pdxt_zy,pdxp_zy,pdyt_zy,pdyp_zy)
pb.plot_efield_beams(pdxt_zy,pdxp_zy,pdyt_zy,pdyp_zy,figno=3)
pb.plot_leakage_beams(l_paper_ideal_zy,rot=[0,0],mn=-0.05,mx=0.05,log=False,figno=4)
plt.savefig('paper_leakage_beams_best.png')

# Let's see the V leakage not vanish - yeah!
l_hera_rich = pb.visibility_leakage_beams(hera_rich['xt'],hera_rich['xp'],hera_rich['yt'],hera_rich['yp'])
pb.plot_leakage_beams(l_hera_rich,rot=[0,90],mn=-3,mx=0,log=True,figno=5)

# And now let's make up the (zero phase) HERA beams at the equator
hdxt_zy = hpt.rotate_healpix_map(np.sqrt(hera_rich['xx']),[0,-90])*pzy['xt']
hdxp_zy = hpt.rotate_healpix_map(np.sqrt(hera_rich['xx']),[0,-90])*pzy['xp']
hdyt_zy = hpt.rotate_healpix_map(np.sqrt(hera_rich['yy']),[0,-90])*pzy['yt']
hdyp_zy = hpt.rotate_healpix_map(np.sqrt(hera_rich['yy']),[0,-90])*pzy['yp']
l_hera_ideal_zy = pb.visibility_leakage_beams(hdxt_zy,hdxp_zy,hdyt_zy,hdyp_zy)
pb.plot_efield_beams(hdxt_zy,hdxp_zy,hdyt_zy,hdyp_zy,figno=6)
pb.plot_leakage_beams(l_hera_ideal_zy,rot=[0,0],mn=-0.05,mx=0.05,log=False,figno=7)
plt.savefig('hera_leakage_beams_best.png')

plt.show()

def fom(leakage):
    qi = (np.power(np.abs(leakage[0,1,:]),2)).sum()
    ui = (np.power(np.abs(leakage[0,2,:]),2)).sum()
    ii = (np.power(np.abs(leakage[0,0,:]),2)).sum()
    f = (qi + ui)/ii
    return f
