"""External libraries used by Where

Currently the package contains the following external library:

+ gpt2w     - Global Pressure and Temperature 2 wet, http://ggosatm.hg.tuwien.ac.at/DELAY/SOURCE/GPT2w/FORTRAN/
+ iers_2010 - Software for IERS Conventions 2010, http://maia.usno.navy.mil/conventions/2010officialinfo.php
+ sofa      - Standards of Fundamental Astronomy, http://www.iausofa.org/


External Fortran libraries are compiled to Python using F2PY. Documentation and examples on how to use F2PY are found
at http://docs.scipy.org/doc/numpy-dev/f2py/.

GPT2W (Global Pressure and Temperature 2 wet)
---------------------------------------------

The GPT2W model is developed by Böhm et.al., and Fortran routines for calculating the model has been made available.


IERS_2010 (Software for IERS Conventions 2010)
----------------------------------------------

Fortran routines corresponding to some of the IERS Conventional models has been made available. See
http://maia.usno.navy.mil/conventions/2010officialinfo.php for more info.


SOFA (Standards of Fundamental Astronomy):
------------------------------------------

The SOFA software library is software maintained by the International Astronomical Union (IAU) for use in astronomical
computing. It is also the reference library used by and referenced by the IERS conventions.

On the webpage http://www.iausofa.org/ Fortran and C versions of the library are available. Using F2PY we have wrapped
the Fortran library and exposed it as a Python module (we wrapped Fortran instead of C as the group is more familiar
with Fortran, and the Fortran library is what has been used in the GEOSAT legacy code). A simple example of its use is
as follows:

    from where.ext import sofa
    from where.data.time import Time

    obsdate = Time(datetime(2015, 8, 7), scale="utc", fmt="datetime")
    earth_rotation_angle = sofa.iau_era00(obsdate.jd, 0)

Documentation and tutorials on how to use the SOFA library can be found at the SOFA webpage,
http://www.iausofa.org/cookbooks.html. Note that the Python wrapping handles return values in a pythonic way. That is,
they are not passed in as arguments to the functions. More precisely, we write

    earth_rotation_angle = sofa.iau_era00(obsdate.jd, 0)

instead of

    sofa.iau_era00(obsdate.jd, 0, earth_rotation_angle)

If you are having problems calling a function properly, do check that the input and output arguments are registered
properly in the `sofa.pyf`-signature file, which is found in the `external/sofa/src`-directory.

IAU regularly updates the SOFA library, both with bugfixes, new data (for instance leap seconds) and improved
models. See the file `external/sofa/UPDATE.txt` for instructions on how to update the Python wrapper.

"""
