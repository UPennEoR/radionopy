# radionopy
Python tools for calculating ionospheric rotation measure behavior for the effects of
the ionosphere for radio telescope measurements.

This is an attempt to make [ionFR](http://sourceforge.net/projects/ionfr) considerably more powerful
and flexible.

The current usage is to edit the executable within rad.py with your latitude and longitude, and the date/time-string. This will generate a folder called "RM_files" with text files of RM healpix maps (hourly) and an npz file containing 24 hours of TEC, RM and RM-uncertainty data. IONEX files required for the date will be obtained automatically via the date/time-string.

`python setup.py develop`

This compiles the [geomagnetic field code](https://www.ngdc.noaa.gov/geomag/WMM/soft.shtml) (within the IGRF folder) 
on the user machine:
`gcc geomag70.c -o geomag70 -lm`


Example Scripts:
-------------
* `ion_altaz.py` -- uses altitudes, azimuths to find RM
* `ion_radec.py` -- uses right ascensions, declinations to find RM

Requirements:
-------------
* [Python](http://www.python.org/) (tested on python 2.7)
* [Numpy](http://scipy.org/)
* [Astropy](http://www.astropy.org/)
* [healpy](http://healpy.readthedocs.org/)

