#!/usr/bin/env python3
"""
Where wrapper for SOFA functions

Description:
------------

This wrapper mainly provides vectorized versions of the SOFA functions.

References:
-----------

.. [1] SOFA Tools for Earth Attitude.
       http://www.iausofa.org/sofa_pn.pdf

.. [2] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
       IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html




"""

# Standard library imports
from functools import lru_cache

# External library imports
import numpy as np

# Where imports
from where.ext import sofa


@lru_cache()
def X_model(time):
    """X coordinates of CIP in GCRS from series based on IAU 2006 precession and IAU 2000A nutation
    """
    return vectorized_xy06(time)[0]


@lru_cache()
def Y_model(time):
    """Y coordinates of CIP in GCRS from series based on IAU 2006 precession and IAU 2000A nutation
    """
    return vectorized_xy06(time)[1]


@lru_cache()
def vectorized_xy06(time):
    """Vectorized version of SOFA xy06-function

    Args:
        time:   lib.time-object

    Returns:
        tuple:  CIP x, y coordinates
    """
    if time.size == 1:
        return sofa.iau_xy06(time.tt.jd_int, time.tt.jd_frac)

    x, y = np.empty(time.shape), np.empty(time.shape)
    for idx, t in enumerate(time.tt):
        x[idx], y[idx] = sofa.iau_xy06(t.jd_int, t.jd_frac)
    return x, y


@lru_cache()
def vectorized_s06(time):
    """Vectorized version of SOFA s06-function

    Args:
        time:   lib.time-object
        x:      CIP x coordinate
        y:      CIP y coordinate

    Returns:
        CIO locator s
    """
    if time.size == 1:
        return sofa.iau_s06(time.tt.jd_int, time.tt.jd_frac, X_model(time), Y_model(time))

    return np.array(
        [sofa.iau_s06(t.jd_int, t.jd_frac, x_i, y_i) for t, x_i, y_i in zip(time.tt, X_model(time), Y_model(time))]
    )


@lru_cache()
def vectorized_era00(time):
    """Vectorized version of SOFA era00-function

    Args:
        time:   lib.time-object

    Returns:
        Earth rotation angle
    """
    if time.size == 1:
        return sofa.iau_era00(time.ut1.jd_int, time.ut1.jd_frac)

    return np.array([sofa.iau_era00(t.jd_int, t.jd_frac) for t in time.ut1])


@lru_cache()
def vectorized_sp00(time):
    """Vectorized version of SOFA sp00-function

    Args:
        time:   lib.time-object

    Returns:
        TIO locator s'
    """
    if time.size == 1:
        return sofa.iau_sp00(time.tt.jd_int, time.tt.jd_frac)

    return np.array([sofa.iau_sp00(t.jd_int, t.jd_frac) for t in time.tt])


@lru_cache()
def vectorized_gmst06(time):
    """Vectorized version of SOFA gmst06-function

    Args:
            Time epochs (see where.data.time for more info)

    Returns:
            Greenwich mean time in radians
    """
    if time.size == 1:
        return sofa.iau_gmst06(time.ut1.jd_int, time.ut1.jd_frac, time.tt.jd_int, time.tt.jd_frac)

    return np.array([sofa.iau_gmst06(t1.jd_int, t1.jd_frac, t2.jd_int, t2.jd_frac) for t1, t2 in zip(time.ut1, time.tt)])


@lru_cache()
def vectorized_gst06(time):
    """Vectorized version of SOFA gmst06-function

    Args:
            Time epochs (see where.data.time for more info)

    Returns:
            Greenwich apparent time in radians
    """

    if time.size == 1:
        return sofa.iau_gst06a(time.ut1.jd_int, time.ut1.jd_frac, time.tt.jd_int, time.tt.jd_frac)

    return np.array([sofa.iau_gst06a(t.ut1.jd_int, t.ut1.jd_frac, t.tt.jd_int, t.tt.jd_frac) for t in time])


@lru_cache()
def vectorized_iau_fal03(time):
    """Vectorized version of SOFA iau_fal03-function

    Args:
           Time epochs (see where.data.time for more info)

    Returns:
           Fundamental argument, mean anomaly of the Moon. In radians.
    """
    julian_centuries = (time.tt.jd - 2_451_545.0) / 36525
    if time.size == 1:
        return sofa.iau_fal03(julian_centuries)
    return np.array([sofa.iau_fal03(t) for t in julian_centuries])


@lru_cache()
def vectorized_iau_falp03(time):
    """Vectorized version of SOFA iau_falp03-function

    Args:
           Time epochs (see where.data.time for more info)

    Returns:
           Fundamental argument, mean anomaly of the Sun. In radians.
    """
    julian_centuries = (time.tt.jd - 2_451_545.0) / 36525
    if time.size == 1:
        return sofa.iau_falp03(julian_centuries)
    return np.array([sofa.iau_falp03(t) for t in julian_centuries])


@lru_cache()
def vectorized_iau_faf03(time):
    """Vectorized version of SOFA iau_faf03-function

    Args:
           Time epochs (see where.data.time for more info)

    Returns:
           Fundamental argument, mean longitude of the Moon minus mean
           longitude of the ascending node. In radians.
    """
    julian_centuries = (time.tt.jd - 2_451_545.0) / 36525
    if time.size == 1:
        return sofa.iau_faf03(julian_centuries)
    return np.array([sofa.iau_faf03(t) for t in julian_centuries])


@lru_cache()
def vectorized_iau_fad03(time):
    """Vectorized version of SOFA iau_fad03-function

    Args:
           Time epochs (see where.data.time for more info)

    Returns:
           Fundamental argument, mean elongation of the Moon from the Sun. In radians.
    """
    julian_centuries = (time.tt.jd - 2_451_545.0) / 36525
    if time.size == 1:
        return sofa.iau_fad03(julian_centuries)
    return np.array([sofa.iau_fad03(t) for t in julian_centuries])


@lru_cache()
def vectorized_iau_faom03(time):
    """Vectorized version of SOFA iau_faom03-function

    Args:
           Time epochs (see where.data.time for more info)

    Returns:
           Fundamental argument, mean longitude of the Moon’s ascending node. In radians.
    """
    julian_centuries = (time.tt.jd - 2_451_545.0) / 36525
    if time.size == 1:
        return sofa.iau_faom03(julian_centuries)
    return np.array([sofa.iau_faom03(t) for t in julian_centuries])
