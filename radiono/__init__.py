'''
radiono

Python tools for calculating ionosphere behavior for the effects of the ionosphere for radio telescope measurements.

Authored by Immanuel Washington, James Aguirre, Saul Kohn, Zac Martinot
2015-2016
'''
import os

#XXX can this be done in a tidier way? Is it required at all?
rad_dir = os.path.dirname(__file__)
root_dir = os.path.abspath(os.path.join(rad_dir, '..'))
rm_dir = os.path.join(root_dir, 'RM_files')
ionex_dir = os.path.join(root_dir, 'TEC')

import utils, ionex_file, interp, physics, rm #<-- I would like this to be the only line of code...
