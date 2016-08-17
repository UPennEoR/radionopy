'''
radiono

Python tools for calculating ionosphere behavior for the effects of the ionosphere for radio telescope measurements.

Authored by Immanuel Washington, James Aguirre, Saul Kohn, Zac Martinot
2015-2016
'''
import os

rad_dir = os.path.dirname(__file__)
root_dir = os.path.abspath(os.path.join(rad_dir, '..'))
rm_dir = os.path.join(root_dir, 'RM_files')
ionex_dir = os.path.join(root_dir, 'TEC')

from rm import RM #XXX I want to remove this call
import utils, ionex_file, interp, physics, rm
