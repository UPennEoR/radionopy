'''
radiono.tests.test_ionex

authors | James Aguirre, Immanuel Washington, Saul Kohn

purpose | Script used to test interp module

Classes
-------
TestIonex | test ionex module
'''
from __future__ import print_function
import os
import unittest
import radiono as rad
from radiono import ionex_file as inx

class TestIonex(unittest.TestCase):
	def setUp(self):
		self.ionex_dir = rad.ionex_dir
		self.year = 2012
		self.month = 4
		self.day = 11
		self.ionex_file = 'CODG1020.12I'
		self.ionex_path = os.path.join(self.ionex_dir, self.ionex_file)
		self.ionex_path_z = ''.join((self.ionex_path, '.z'))

	def tearDown(self):
		pass

	def test_IONEX_file(self):
		ionex_path = inx.IONEX_file_needed(self.year, self.month, self.day, self.ionex_dir)
		self.assertEqual(ionex_path, self.ionex_path, msg='IONEX file downloaded incorrect')

		ionex_path_z = inx.get_IONEX_file(self.year, self.month, self.day, self.ionex_path, self.ionex_dir)
		self.assertEqual(ionex_path_z, self.ionex_path, msg='IONEX file downloaded incorrect')

	def test_gen_IONEX_list(self):
		pass

	def test_read_IONEX_TEC(self):
		pass

	def test_IONEX_data(self):
		pass

if __name__ == '__main__':
	unittest.main()