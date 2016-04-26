from setuptools import find_packages
import os
import subprocess

__version__ = '1.0.0-dev'

setup_args = {
    'name': 'radionopy',
    'author': 'Immanuel Washington',
    'author_email': 'immwa at sas.upenn.edu',
    'description': 'package for ionosphere RM',
    'url': 'https://github.com/jaguirre/radionopy.git',
    'license': '?',
    'package_dir' : {'radionopy': ''},
    'packages' : find_packages(),
    'version': __version__,
}

base_path = os.path.expanduser(os.getcwd())
rad_path = os.path.join(base_path, 'radiono')

if __name__ == '__main__':
    try:
        from setuptools import setup
    except:
        from distutils.core import setup
    try:
        apply(setup, (), setup_args)
    except:
        setup(**setup_args)

    #Compiles geomag C script
    script_program = 'gcc'
    script_data = os.path.join(rad_path, 'IGRF/geomag70_linux/geomag70.c')
    script_name = os.path.join(rad_path, 'IGRF/geomag70_linux/geomag70')
    script_option = '-o'
    script_fix = '-lm'
    subprocess.call([script_program, script_data, script_option, script_name, script_fix])
