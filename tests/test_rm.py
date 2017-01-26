'''
radiono.tests.test_rm

authors | Immanuel Washington, Saul Kohn

purpose | Script used to test IonoMap class in rm.py
'''
from __future__ import print_function
import unittest, os, shutil, random
from radiono import rm
from radiono import utils as ut
import healpy as hp, numpy as np, matplotlib.pyplot as plt

paperLocStrings = ('30d43m17.5ss','21d25m41.9se')
lofarLocStrings = ('53d00m00.00sn','6d00m00.00se')
testIonexDir = './localIonex'
testRmDir = './localRM'
testHeight = 1000
testTime='2011-04-11'
test_yr,test_month,test_day = map(int,testTime.split('-'))
#Cas A
testRA = np.radians((23.+(23./60.)+(27.9/3600.))*15.)
testDec= np.radians(58.+(48./60.)+(42.4/3600.))

testRAs,testDecs = ut.nsideToRaDec(16)
npix = hp.nside2npix(16)

class TestIonoMap(unittest.TestCase):
    def setUp(self):
        self.rm_map = rm.IonoMap(paperLocStrings[0], paperLocStrings[1], [testTime], height=testHeight, ionex_dir=testIonexDir, rm_dir=testRmDir)
    def test_map_properties(self):
        self.assertEqual(self.rm_map.lat.value, -30.721527777777776)
        self.assertEqual(self.rm_map.lon.value, 21.428305555555557)
        self.assertEqual(self.rm_map.height, testHeight)
        self.assertEqual(np.around(self.rm_map.location.latitude.value,5), np.around(-30.721527777777776,5))
        self.assertEqual(np.around(self.rm_map.location.longitude.value,5), np.around(21.428305555555557,5))
        self.assertEqual(np.around(self.rm_map.location.height.value,5), np.around(testHeight,5))

        self.assertEqual(self.rm_map.rm_dir, testRmDir)
        self.assertEqual(self.rm_map.ionex_dir, testIonexDir)
    """
    def test_make_rm_dir(self):
        self.rm_map.make_rm_dir(testTime)
        assert(os.path.exists(testRmDir+'/%s'%testTime))
    """
    def test_ionex_data(self):
        tec_hp,rms_hp,ion_height = self.rm_map.ionex_data(test_yr,test_month,test_day,ionex_dir=testIonexDir)
        assert(tec_hp.shape == rms_hp.shape)
        assert(tec_hp.shape[0]==24)
        hp.npix2nside(tec_hp.shape[1])
        assert(ion_height > self.rm_map.location.height.value)

    def test_altaz(self):
        rm_altaz_map = self.rm_map.altaz()

    def test_radec_single(self):
        self.lofar_map = rm.IonoMap(lofarLocStrings[0],lofarLocStrings[1], [testTime], height=0)
        self.lofar_map.get_radec_RM([testRA],[testDec])
        assert(self.lofar_map.RMs.shape == (1,24,1)) #npix = 1
        assert(self.lofar_map.dRMs.shape == (1,24,1))
        #Plot the comparison plot for Sotomayor Beltran et al. 2013 Fig.4b
        plt.errorbar(range(24), np.abs(self.lofar_map.RMs[0,:,0]), yerr=self.lofar_map.dRMs[0,:,0], fmt='ro', ecolor='r')
        plt.xlim(0,23)
        plt.ylim(0,2.5)
        plt.xlabel(r'UT Time [Hours]')
        plt.ylabel(r'$|\phi_{\rm ion}|$ [rad m$^{-2}$]')
        plt.savefig('testFig_compare_SB_Fig4b.png')
        plt.suptitle(r'%s'%('/'.join(map(str,[test_day,test_month,test_yr]))))
        plt.close()

    def test_radec_multip(self):
        self.rm_map.get_radec_RM(testRAs,testDecs)
        assert(self.rm_map.RMs.shape == (1,24,npix))
        assert(self.rm_map.dRMs.shape == (1,24,npix))
        hp.orthview(self.rm_map.RMs[0,0,:])
        plt.close()

    def tearDown(self):
        if os.path.exists(testRmDir): shutil.rmtree(testRmDir)
        if os.path.exists(testIonexDir): shutil.rmtree(testIonexDir)
        return None


if __name__ == '__main__':
    unittest.main()
