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
from where.lib import constant
from where.ext import sofa
from where.lib.unit import unit
from where.lib import rotation
from where import apriori


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
    return rotation.R3(-E(time)) @ rotation.R2(-d(time)) @ rotation.R3(E(time)) @ rotation.R3(s(time))

@cache.function
def s(time):
    """CIO locator
    """
    return vectorized_s06(time)

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
def X(time):
    """Total X coordinate of CIP

    Returns:
        model + celestial pole offset
    """
    eop = apriori.get("eop", time=time)
    return X_model(time) + eop.dx * unit.arcsec2rad

@cache.function
def Y(time):
    """Total Y coordinate of CIP

    Returns:
        model + celestial pole offset
    """
    eop = apriori.get("eop", time=time)
    return Y_model(time) + eop.dy * unit.arcsec2rad

@cache.function
def Z(time):
    """Total Z coordinate of CIP

    Derived from X and Y with assumption of unit sphere with

    Returns:
        model + celestial pole offset
    """
    return np.sqrt(1 - (X(time)**2 + (Y(time)**2)))

@cache.function
def d(time):
    """Polar angle of CIP from positive z-axis

    Spherical coordinate with radius = 1
    """
    return np.arccos(Z(time))

@cache.function
def E(time):
    """Azimuth angle of CIP in xy-plane

    Spherical coordinate with radius = 1
    """
    return np.arctan2(Y(time), X(time))

@cache.function
def ERA(time):
    """Earth rotation angle
    """
    return vectorized_era00(time)

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

    return rotation.R3(-ERA(time))


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

    return -rotation.dR3(-ERA(time)) * constant.omega

@cache.function
def xp(time):
    """X coordinate of the CIP in ITRS
    """
    eop = apriori.get("eop", time=time)
    return eop.x * unit.arcsec2rad

@cache.function
def yp(time):
    """Y coordinate of the CIP in ITRS
    """
    eop = apriori.get("eop", time=time)
    return eop.y * unit.arcsec2rad

@cache.function
def s_prime(time):
    """ TIO locator
    """
    return vectorized_sp00(time)


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
    return rotation.R3(-s_prime(time)) @ rotation.R2(xp(time)) @ rotation.R1(yp(time))

@cache.function
def dW_dxp(time):
    """Derivative of transformation matrix for the polar motion with regards to the CIP (Celestial Intermediate Pole)
    in TRF along the Greenwich meridian x_p.

    This is done according to equations (2.31) in Teke :cite:`teke2011` which is the analytical partial derivative of
    equation (5.3) in :cite:'iers2010'.
    """
    return rotation.R3(-s_prime(time)) @ rotation.dR2(xp(time)) @ rotation.R1(yp(time))


@cache.function
def dW_dyp(time):
    """Derivative of transformation matrix for the polar motion with regards to the CIP in ITRS.

    Analytical partial derivative of equation (5.3) in :cite:'iers2010'.
    """
    return rotation.R3(-s_prime(time)) @ rotation.R2(xp(time)) @ rotation.dR1(yp(time))

@cache.function
def dE_dX(time):
    """Derivative of azimuth angle of CIP in GCRS with regards to the X coordinate of CIP in GCRS
    """
    return (-Y(time) / (X(time) ** 2 + Y(time) ** 2))[:, None, None]

@cache.function
def ds_dX(time):
    """Derivative of CIO locator with regards the X coordinate of CIP in GCRS
    """
    return (-Y(time) / 2)[:, None, None]

@cache.function
def dd_dX(time):
    """Derivative of polar angle of CIP in GCRS with regards to the X coordinate of CIP in GCRS
    """
    return (X(time) / (Z(time) * np.sqrt(X(time) ** 2 + Y(time) ** 2)))[:, None, None]

@cache.function
def dE_dY(time):
    """Derivative of azimuth angle of CIP in GCRS with regards to the Y coordinate of CIP in GCRS
    """
    return (X(time) / (X(time) ** 2 + Y(time) ** 2))[:, None, None]

@cache.function
def ds_dY(time):
    """Derivative of CIO locator with regards the X coordinate of CIP in GCRS
    """
    return (-X(time) / 2)[:, None, None]

@cache.function
def dd_dY(time):
    """Derivative of polar angle of CIP in GCRS with regards to the Y coordinate of CIP in GCRS
    """
    return (Y(time) / (Z(time) * np.sqrt(X(time) ** 2 + Y(time) ** 2)))[:, None, None]

@cache.function
def dQ_dX(time):
    """Derivative of transformation matrix for nutation/presession with regards to the X coordinate of CIP in GCRS
    """
    # Rotation matrices
    R3_E = rotation.R3(E(time))
    R3_s = rotation.R3(s(time))
    R2_md = rotation.R2(-d(time))
    R3_mE = rotation.R3(-E(time))
    dR3_s = rotation.dR3(s(time))
    dR3_E = rotation.dR3(E(time))
    dR3_mE = rotation.dR3(-E(time))
    dR2_md = rotation.dR2(-d(time))

    return (
          dR3_mE @  R2_md @  R3_E @  R3_s * (-dE_dX(time))
        +  R3_mE @ dR2_md @  R3_E @  R3_s * (-dd_dX(time))
        +  R3_mE @  R2_md @ dR3_E @  R3_s * ( dE_dX(time))
        +  R3_mE @  R2_md @  R3_E @ dR3_s * ( ds_dX(time))
    )

@cache.function
def dQ_dY(time):
    """Derivative of transformation matrix for nutation/presession with regards to the Y coordinate of CIP in GCRS
    """
    # Rotation matrices
    R3_E = rotation.R3(E(time))
    R3_s = rotation.R3(s(time))
    R2_md = rotation.R2(-d(time))
    R3_mE = rotation.R3(-E(time))
    dR3_s = rotation.dR3(s(time))
    dR3_E = rotation.dR3(E(time))
    dR3_mE = rotation.dR3(-E(time))
    dR2_md = rotation.dR2(-d(time))

    return (
          dR3_mE @  R2_md @  R3_E @  R3_s * (-dE_dY(time))
        +  R3_mE @ dR2_md @  R3_E @  R3_s * (-dd_dY(time))
        +  R3_mE @  R2_md @ dR3_E @  R3_s * ( dE_dY(time))
        +  R3_mE @  R2_md @  R3_E @ dR3_s * ( ds_dY(time))
    )

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

    return np.array([sofa.iau_s06(t.jd1, t.jd2, x_i, y_i) for t, x_i, y_i in zip(time.tt, X_model(time), Y_model(time))])


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
