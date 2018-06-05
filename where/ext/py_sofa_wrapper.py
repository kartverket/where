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



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

# Standard library imports

# External library imports
import numpy as np

# Where imports
from where.lib import cache
from where.lib import constant
from where.ext import sofa
from where.lib.unit import unit
from where.lib import rotation


@cache.function
def Q(time):
    """Transformation matrix for the celestial motion of the CIP

    The transformation matrix is described in the IERS Conventions [2], section 5.4.4. The implementation is based on
    section 5.6 "IAU 2006/2000A, CIO based, using X, Y series" in the SOFA Tools for Earth Attitude [1].

    Args:
        time:  A lib.time Time-object

    Returns:
        numpy.array: An array of 3x3 transformation matrices
    """
    # Read EOP values
    from where import apriori

    eop = apriori.get("eop", time=time)

    # CIP and CIO
    x, y = vectorized_xy06(time)
    s = vectorized_s06(time, x, y)

    # Add Celestial Intermediate Pole corrections
    x += eop.dx * unit.arcsec2rad
    y += eop.dy * unit.arcsec2rad

    # Celestial pole motion
    if time.isscalar:
        return vectorized_c2ixys(x, y, s).T
    else:
        return vectorized_c2ixys(x, y, s).transpose(0, 2, 1)


@cache.function
def R(time):
    """Transformation matrix for the Earth rotation

    The transformation matrix is described in the IERS Conventions [2]_, section 5.4.2. The implementation is based on
    section 5.6 "IAU 2006/2000A, CIO based, using X, Y series" in the SOFA Tools for Earth Attitude [1]_.

    According to equation (5.5) in the IERS Conventions [2]_, `The CIO based transformation matrix arising from the
    rotation of the Earth around the axis of the CIP can be expressed as`

    .. math::

        R(t) = R_3(-ERA)

    Args:
        time:  A lib.time Time-object

    Returns:
        numpy.array: An array of 3x3 transformation matrices

    """
    era = vectorized_era00(time)
    return rotation.R3(-era)


@cache.function
def dR_dut1(time):
    """Derivative of transformation matrix for the Earth rotation with respect to time (UT1??)

    The transformation matrix is described in the IERS Conventions [2]_, section 5.4.2. The implementation is based on
    section 5.6 "IAU 2006/2000A, CIO based, using X, Y series" in the SOFA Tools for Earth Attitude [1]_.

    Unit of return value: radians / seconds ??

    Args:
        time:  A lib.time Time-object

    Returns:
        numpy.array: An array of 3x3 transformation matrices

    """
    era = vectorized_era00(time)
    return -rotation.dR3(-era) * constant.omega


@cache.function
def W(time):
    """Transformation matrix for the polar motion

    The transformation matrix is described in the IERS Conventions [2], section 5.4.1. The implementation is based on
    section 5.6 "IAU 2006/2000A, CIO based, using X, Y series" in the SOFA Tools for Earth Attitude [1].

    Args:
        time:  A lib.time Time-object

    Returns:
        numpy.array: An array of 3x3 transformation matrices
    """
    # Read EOP values
    from where import apriori

    eop = apriori.get("eop", time=time)

    # Polar motion
    xp = eop.x * unit.arcsec2rad
    yp = eop.y * unit.arcsec2rad
    sp = vectorized_sp00(time)

    if time.isscalar:
        return vectorized_pom00(xp, yp, sp).T
    else:
        return vectorized_pom00(xp, yp, sp).transpose(0, 2, 1)


@cache.function
def dW_dxp(time):
    """Derivative of transformation matrix for the polar motion with regards to the CIP (Celestial Intermediate Pole)
    in TRF along the Greenwich meridian x_p.

    This is done according to equations (2.31) in Teke :cite:`teke2011` which is the analytical partial derivative of
    equation (5.3) in :cite:'iers2010'.
    """
    from where import apriori

    eop = apriori.get("eop", time=time)

    xp = eop.x * unit.arcsec2rad
    yp = eop.y * unit.arcsec2rad

    sp = vectorized_sp00(time)

    R3_sp = rotation.R3(sp)
    R1_yp = rotation.R1(yp)
    dR2_xp = rotation.dR2(xp)

    return R3_sp @ dR2_xp @ R1_yp


@cache.function
def dW_dyp(time):
    """Derivative of transformation matrix for the polar motion with regards to the CIP (Celestial Intermediate Pole)
    in TRF along 270degrees longitude y_p.

    This is done according to equations (2.33) in Teke :cite:`teke2011` which is the analytical partial derivative of
    equation (5.3) in :cite:'iers2010'.
    """
    from where import apriori

    eop = apriori.get("eop", time=time)

    xp = eop.x * unit.arcsec2rad
    yp = eop.y * unit.arcsec2rad

    sp = vectorized_sp00(time)

    R3_sp = rotation.R3(sp)
    R2_xp = rotation.R2(xp)
    dR1_yp = rotation.dR1(yp)

    return R3_sp @ R2_xp @ dR1_yp


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


def vectorized_s06(time, x, y):
    """Vectorized version of SOFA s06-function

    Args:
        time:   lib.time-object
        x:      CIP x coordinate
        y:      CIP y coordinate

    Returns:
        CIO locator s
    """
    if time.isscalar:
        return sofa.iau_s06(time.tt.jd1, time.tt.jd2, x, y)

    return np.array([sofa.iau_s06(t.jd1, t.jd2, x_i, y_i) for t, x_i, y_i in zip(time.tt, x, y)])


def vectorized_c2ixys(x, y, s):
    """Vectorized version of SOFA c2ixys-function

    Args:
        x:      CIP x coordinate.
        y:      CIP y coordinate.
        s:      CIO locator.

    Returns:
        np.array:  Celestial to intermediate matrices
    """
    if isinstance(s, (float, int)):
        return sofa.iau_c2ixys(x, y, s)
    return np.array([sofa.iau_c2ixys(x_i, y_i, s_i) for x_i, y_i, s_i in zip(x, y, s)])


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


def vectorized_pom00(xp, yp, sp):
    """Vectorized version of SOFA pom00-function

    Args:
        xp:      CIP x polar coordinate.
        yp:      CIP y polar coordinate.
        sp:      TIO locator.

    Returns:
        np.array:  Polar motion matrix
    """
    if isinstance(sp, (float, int)):
        return sofa.iau_pom00(xp, yp, sp)

    return np.array([sofa.iau_pom00(x_i, y_i, s_i) for x_i, y_i, s_i in zip(xp, yp, sp)])


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
