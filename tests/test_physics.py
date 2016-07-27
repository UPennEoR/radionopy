'''
radiono.tests.test_physics

authors | James Aguirre, Immanuel Washington, Saul Kohn

purpose | Script used to test physics module

Classes
-------
TestPhysics | test physics module
'''
from __future__ import print_function
import unittest
import radiono as rad
from radiono import physics as phys

class TestPhysics(unittest.TestCase):
	def setUp(self):
	    self.nside = 16
	    self.npix = hp.nside2npix(self.nside)
	    self.ipix = np.arange(self.npix)
	    self.theta, self.phi = hp.pix2ang(self.nside, self.ipix)

	    self.alt_src = 90. - np.degrees(np.array(self.theta))
	    self.az_src = np.degrees(np.array(self.phi))
		self.lat_str = '30d43m17.5ss'
	    self.lon_str = '21d25m41.9se'
	    self.dates = (('2012-11', range(18, 31)),
	        	      ('2012-12', range(1, 32)),
	            	  ('2013-01', range(1, 32)),
	           		  ('2013-02', range(1, 29)),
	           		  ('2013-03', range(1, 32)))

	    self.ra_strs = ('16h50m04.0s',)
	    self.dec_strs = ('+79d11m25.0s',)
	    
	    self.time_part = 'T00:00:00'
	    self.date_strs = ('-'.join((ym, rad.std_hour(day, verbose=False))) for ym, days in self.dates for day in days)
	    self.time_strs = (''.join((date_str, self.time_part)) for date_str in self.date_strs)

	    self.height = 0

	def tearDown(self):
		pass

	def test_B_IGRF(self):
		pass

	def test_punct_ion_offset(self):
		pass

	def test_get_coords(self):
	    #lat_obs = Latitude(Angle(lat_str[:-1]))
	    #lon_obs = Longitude(Angle(lon_str[:-1]))
		#coord_lat, coord_lon = phys.get_coords(self.lat_str, self.lon_str, lat_obs, lon_obs, off_lat, off_lon)
		pass

	def test_ipp(self):
		#ipp = phys.ipp(self.lat_str, self.lon_str, self.az_src, 90. - self.alt_src, self.height)
		pass

	def test_rotate_healpix_map(self):
		pass

if __name__ == '__main__':
	unittest.main()