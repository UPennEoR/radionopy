'''
radiono
'''
import os
from astropy import constants as c

rad_dir = os.path.dirname(__file__)
root_dir = os.path.abspath(os.path.join(rad_dir, '..'))
rm_dir = os.path.join(root_dir, 'RM_files')
ionex_dir = os.path.join(root_dir, 'TEC')
TECU = 1e16
TEC2m2 = 0.1 * TECU
earth_radius = c.R_earth.value #6371000.0 # in meters
tesla_to_gauss = 1e4

from rm import RM #XXX I want to remove this call
import utils, ionex_file, interp, physics, rm
