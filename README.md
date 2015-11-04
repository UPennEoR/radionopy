# radionopy
Python tools for calculating ionosphere behavior for the effects of
the ionosphere for radio telescope measurements.

This is an attempt to make ionFR
(http://sourceforge.net/projects/ionfr/) considerably more powerful
and flexible.

The general flow goes like this, and I have tried indicate where in
the code the various bits are currently implemented.

- determine which IONEX file you need for a given day of interest.
  Currently implemented in getIONEXfile.py

- generate the TEC map from the IONEX file.  This is done in readTEC.py.

  readIonexTEC will correctly return 13 maps in lat / lon coordinates, but
  does not implement the time interpolation (which is slightly non-trivial).
  I did not use the the spatial interpolation scheme; my idea was to grid onto
  a HealPix map.
  (This is based on the original function calcTEC in ionFR, a hacked version
  of which is present in readTEC.py)

- calculate the "pierce points" for every direction of interest

- calculate the B-field at the pierce points

  This depends on using the IGRF code, which is a standalone C
  executable which requires the generation of an input file and
  reading an output file.  For now, I think we find a sensible way to
  wrap this, but should look if it is the limiting step for
  calculating a large number of points.

- calculate vertical TEC at the pierce points.  This uses 

 