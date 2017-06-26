# radionopy
[![Build Status](https://travis-ci.org/UPennEoR/radionopy.svg?branch=master)](https://travis-ci.org/UPennEoR/radionopy)
[![Coverage Status](https://coveralls.io/repos/github/UPennEoR/radionopy/badge.svg?branch=master)](https://coveralls.io/github/UPennEoR/radionopy?branch=master)

Python tools for calculating ionospheric rotation measure behavior for the effects of
the ionosphere for radio telescope measurements.

This is an attempt to make [ionFR](http://sourceforge.net/projects/ionfr) considerably more powerful
and flexible.

Installation:
-------------
`python setup.py develop`

This compiles the [geomagnetic field code](https://www.ngdc.noaa.gov/geomag/WMM/soft.shtml) (within the IGRF folder) 
on the user machine:
`gcc geomag70.c -o geomag70 -lm`


Example Scripts:
-------------
* `docs/radionopy_Tutorial.ipynb` -- tutorial in a Jupyter notebook.
* `scripts/probe_parameter_space.py` -- creating and interacting with an IonoMap object containing months-worth of data.

Requirements:
-------------
* [Python](http://www.python.org/) (tested on python 2.7)
* [Numpy](http://scipy.org/)
* [Astropy](http://www.astropy.org/)
* [healpy](http://healpy.readthedocs.org/)
* [pyEphem](https://pypi.python.org/pypi/ephem/)
* optional: [Basemap](http://matplotlib.org/basemap/users/index.html)
