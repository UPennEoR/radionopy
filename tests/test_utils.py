from radiono import utils as ut
import unittest, os, shutil

testDatedash = '2011-01-05'
testDateslash = '2011/01/05'

class TestUtils(unittest.TestCase):
    def test_stdHour(self):
        hr1 = ut.std_hour(1)
        self.assertEqual(hr1,'01') #should be "01"
        hr10 = ut.std_hour(10,verbose=True)
    def test_eph2ionDate(self):
        ionDate = ut.eph2ionDate(testDatedash)
        self.assertEqual(ionDate,testDatedash)
    def test_ion2ephDate(self):
        ephDate = ut.ion2ephDate(testDateslash)
        self.assertEqual(ephDate,testDatedash)
    def test_ephemPAPER(self):
        site = ut.ephemPAPER(date=testDateslash)
        self.assertEqual(site.elevation,1000.0)
        self.assertEqual(site.date,40546.5)
        self.assertEqual(site.lon,0.37399448506783706)
        self.assertEqual(site.lat,-0.5361918109651189)
    def test_nextTransit(self):
        pass
    def test_parseTransitBasic(self):
        pass 
