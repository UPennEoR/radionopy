# radionopy
Python tools for calculating ionosphere behavior for the effects of
the ionosphere for radio telescope measurements.

This is an attempt to make ionFR
(http://sourceforge.net/projects/ionfr/) considerably more powerful
and flexible.

The current usage is to edit the executable within rad.py with your latitude and longitude, and the date/time-string. This will generate a folder called "RM files" with text files of RM healpix maps (hourly) and an npz file containing 24 hours of TEC, RM and uRM data. IONEX files required for the date will be obtained automatically via the date/time-string.

In order to execute this successfully, one must compile the geomagnetic field C code (within the IGRF folder) on their own machine:

`gcc geomag70.c -o geomag70`

Main Scripts:
-------------
* `radiono.py` -- main module
* `ion_altaz.py` -- uses altitudes, azimuths to find RM
* `ion_radec.py` -- uses right ascensions, declinations to find RM

Requirements:
-------------
* `Python <http://www.python.org/>`
* `Numpy <http://scipy.org/>`
* `Astropy <http://www.astropy.org/>`
* `healpy <http://healpy.readthedocs.org/>`

Optional:
---------
* `PyGSM <https://github.com/telegraphic/PyGSM>`
* `pyephem <https://github.com/brandon-rhodes/pyephem/>`

