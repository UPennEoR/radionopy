'''
radiono.tests.test_rm

authors | James Aguirre, Immanuel Washington, Saul Kohn

purpose | Script used to test RM class

Classes
-------
TestRM | test RM class module
'''
from __future__ import print_function
import unittest, os, shutil
from radiono import rm
import healpy as hp, numpy as np
from matplotlib import pyplot as plt

paperLocStrings = ('30d43m17.5ss','21d25m41.9se')
testIonexDir = './localIonex'
testRmDir = './localRM'
testHeight = 1000
testNside=16
testTime='2010-01-21T00:00:00'
testRA = ['0.0']
testDec= ['-30.7']

class TestRM(unittest.TestCase):
    def setUp(self):
        self.map = rm.RM(paperLocStrings[0], paperLocStrings[1], [testTime], height=testHeight,\
        nside=testNside, ionex_dir=testIonexDir, rm_dir=testRmDir)
    def test_map_properties(self):
        self.assertEqual(self.map.lat.value, -30.721527777777776)
        self.assertEqual(self.map.lon.value, 21.428305555555557)
        self.assertEqual(self.map.height, testHeight)
        
        self.assertEqual(np.around(self.map.location.latitude.value,5), np.around(-30.721527777777776,5))
        self.assertEqual(np.around(self.map.location.longitude.value,5), np.around(21.428305555555557,5))
        self.assertEqual(np.around(self.map.location.height.value,5), np.around(testHeight,5))    
    
        self.assertEqual(self.map.rm_dir, testRmDir)
        self.assertEqual(self.map.ionex_dir, testIonexDir)
    
    def test_make_rm_dir(self):
        self.map.make_rm_dir(testTime,verbose=True)
        assert(os.path.exists(testRmDir+'/2010-01-21'))
    
    def test_ionex_data(self):
        tec_hp,rms_hp,ion_height = self.map.ionex_data(2010,01,21,ionex_dir=testIonexDir,verbose=True)
        assert(tec_hp.shape == rms_hp.shape)
        assert(tec_hp.shape[0]==24)
        hp.npix2nside(tec_hp.shape[1])
        assert(ion_height > self.map.location.height.value)
        
    def tearDown(self):
        if os.path.exists(testRmDir): shutil.rmtree(testRmDir)
        if os.path.exists(testIonexDir): shutil.rmtree(testIonexDir)
        


if __name__ == '__main__':
    unittest.main()
