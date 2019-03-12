#!/usr/bin/env python3
"""
Where wrapper for SOFA functions

Description:
------------

This wrapper mainly provides vectorized versions of the SOFA functions, as well as cached versions of the most
important matrices like Q (the transformation matrix for the celestial motion of the CIP), R (the transformation matrix
for Earth rotation) and W (the transformation matrix for polar motion). See the IERS conventions [2]_ for more details,
in particular section 5.

References:
-----------

.. [1] SOFA Tools for Earth Attitude.
       http://www.iausofa.org/sofa_pn.pdf

.. [2] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
       IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html




"""

# Standard library imports

# External library imports
import numpy as np

# Where imports
from where.lib import cache
from where.ext import sofa


@cache.function
def X_model(time):
    """X coordinates of CIP in GCRS from series based on IAU 2006 precession and IAU 2000A nutation
    """
    return vectorized_xy06(time)[0]


@cache.function
def Y_model(time):
    """Y coordinates of CIP in GCRS from series based on IAU 2006 precession and IAU 2000A nutation
    """
    return vectorized_xy06(time)[1]


@cache.function
def vectorized_xy06(time):
    """Vectorized version of SOFA xy06-function

    Args:
        time:   lib.time-object

    Returns:
        tuple:  CIP x, y coordinates
    """
    if time.isscalar:
        return sofa.iau_xy06(time.tt.jd1, time.tt.jd2)

    x, y = np.empty(time.shape), np.empty(time.shape)
    for idx, t in enumerate(time.tt):
        x[idx], y[idx] = sofa.iau_xy06(t.jd1, t.jd2)
    return x, y


@cache.function
def vectorized_s06(time):
    """Vectorized version of SOFA s06-function

    Args:
        time:   lib.time-object
        x:      CIP x coordinate
        y:      CIP y coordinate

    Returns:
        CIO locator s
    """
    if time.isscalar:
        return sofa.iau_s06(time.tt.jd1, time.tt.jd2, X_model(time), Y_model(time))

    return np.array(
        [sofa.iau_s06(t.jd1, t.jd2, x_i, y_i) for t, x_i, y_i in zip(time.tt, X_model(time), Y_model(time))]
    )


@cache.function
def vectorized_era00(time):
    """Vectorized version of SOFA era00-function

    Args:
        time:   lib.time-object

    Returns:
        Earth rotation angle
    """
    if time.isscalar:
        return sofa.iau_era00(time.ut1.jd1, time.ut1.jd2)

    return np.array([sofa.iau_era00(t.jd1, t.jd2) for t in time.ut1])


@cache.function
def vectorized_sp00(time):
    """Vectorized version of SOFA sp00-function

    Args:
        time:   lib.time-object

    Returns:
        TIO locator s'
    """
    if time.isscalar:
        return sofa.iau_sp00(time.tt.jd1, time.tt.jd2)

    return np.array([sofa.iau_sp00(t.jd1, t.jd2) for t in time.tt])


def vectorized_llh(pos, ref_ellipsoid=2):
    """Vectorized version of SOFA gc2gd-function. 

    Converts xyz coordinates to latitude, longitude and height

    @todo ref_ellipsoid
    Args:
        pos:        xyz coordinates
    Returns:
        np.array:   llh coordinates 
    """
    if np.array(pos).ndim == 1:
        lon, lat, h, _ = sofa.iau_gc2gd(ref_ellipsoid, pos)
        return np.array([lat, lon, h])
    else:
        llh = np.empty(pos.shape)
        for i, xyz in enumerate(pos):
            lon, lat, h, _ = sofa.iau_gc2gd(ref_ellipsoid, xyz)
            llh[i, :] = (lat, lon, h)
        return llh
