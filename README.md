# radionopy
Python tools for calculating ionosphere behavior for the effects of
the ionosphere for radio telescope measurements.

This is an attempt to make ionFR
(http://sourceforge.net/projects/ionfr/) considerably more powerful
and flexible.

The general flow goes like this, based on ionFRM.py, and I have tried
indicate where in the code the various bits are currently implemented.

- determine which IONEX file you need for a given day of interest.
  Currently implemented in getIONEXfile.py

- calculate the "pierce points" for every direction of interest.
  Defined in the function PuncIonOffset in ippcoor.py.  Note that this
  uses math instead of numpy, so is not vectorized.

- generate the TEC map from the IONEX file.  This is done in readTEC.py.

  readIonexTEC will correctly return 13 maps in lat / lon coordinates,
  but does not implement the time interpolation (which is slightly
  non-trivial).  I did not use the the spatial interpolation scheme;
  my idea was to grid onto a HealPix map.  I would like to keep the
  whole maps around, unlike ionFR, which was only interested in
  particular positions.
  
  This is based on the original function calcTEC in IONEX/teccalc.py
  in ionFR, a hacked version of which is present in readTEC.py.  It
  appears that RMS errors on TEC were obtained with the almost
  identical tecrmscalc ... clearly these should be merged into one
  function.

- there is an important division by cos(ZenPunct) which is not
  embedded in a fuction

- calculate the B-field at the pierce points

  This depends on using the IGRF code, which is a standalone C
  executable which requires the generation of an input file and
  reading an output file.  For now, I think we find a sensible way to
  wrap this, but should look if it is the limiting step for
  calculating a large number of points.

  f = open(''+str(path)+'ionFR/IGRF/geomag70_linux/input.txt', 'w')

- dot the B field into the direction

 