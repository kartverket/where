"""Library for geodetic rotation matrices

Description:
------------

Creates rotation matrices for rotation around the axes of a right handed Cartesian coordinate system and
their derivatives.

Contains the rotation matrices (and angles) for the transition from a terrestrial reference system to a 
geocentric celestial reference system according to the IERS 2010 Conventions.
"""
# Standard library imports
from functools import lru_cache

# Include all rotation matrices defined in Midgard
from midgard.math.constant import constant
from midgard.math.rotation import *  # noqa
from midgard.math.unit import Unit

# Where imports
from where import apriori
from where.ext import sofa_wrapper as sofa

#
# Transformations
#


@lru_cache()
def gcrs2trs(time):
    """Transformation from space fixed to earth fixed coordinate system

    According to IERS 2010 Conventions
    """
    if time.size == 1:
        return (Q(time) @ R(time) @ W(time)).transpose()
    else:
        return (Q(time) @ R(time) @ W(time)).transpose(0, 2, 1)


@lru_cache()
def trs2gcrs(time):
    """Transformation from earth fixed to space fixed coordinate system

    According to IERS 2010 Conventions
    """
    return Q(time) @ R(time) @ W(time)


@lru_cache()
def dtrs2gcrs_dt(time):
    """Derivative of transformation from earth fixed to space fixed coordinate system with regards to time

    dQ/dt and dW/dt is approximated to zero since these matrices change slowly compared to dR/dt 
    """
    return Q(time) @ dR_dut1(time) @ W(time)


@lru_cache()
def dgcrs2trs_dt(time):
    """Transformation from space fixed to earth fixed coordinate system

    dQ/dt and dW/dt is approximated to zero since these matrices change slowly compared to dR/dt
    """
    if time.size == 1:
        return (Q(time) @ dR_dut1(time) @ W(time)).transpose()
    else:
        return (Q(time) @ dR_dut1(time) @ W(time)).transpose(0, 2, 1)


# @lru_cache() # TODO sat_pos unhashable
def yaw2trs(sat_pos, time):
    """Transformation matrix from yaw-steering reference system to ITRS."""
    return trs2yaw(sat_pos, time).transpose(0, 2, 1)


# @lru_cache() # TODO sat_pos unhashable
def trs2yaw(sat_pos, time):
    """Transformation matrix from ITRS to yaw-steering reference system

    The yaw-steering reference system given with x-axis lying in the Earth-Satellite-Sun plane, y-axis as the
    normal vector of the Earth-Satellite-Sun plane and the z-axis pointing to the Earth's center.
    """
    eph = apriori.get("ephemerides", time=time.tdb)
    z_unit = -sat_pos.trs.pos.unit_vector  # unit vector of z-axis
    sat_sun = eph.pos_itrs("sun") - sat_pos.trs.pos  # vector pointing from satellite position to Sun
    y = np.cross(z_unit, sat_sun)
    y_unit = y / np.linalg.norm(y, axis=1)[:, None]  # unit vector of y-axis
    x = np.cross(y_unit, z_unit)
    x_unit = x / np.linalg.norm(x, axis=1)[:, None]  # unit vector of z-axis

    return np.stack((x_unit, y_unit, z_unit), axis=1)


#
# Rotation matrices
#
@lru_cache()
def Q(time):
    """Transformation matrix for the celestial motion of the CIP

    The transformation matrix is described in the IERS Conventions [2], section 5.4.4. The implementation is based on
    section 5.6 "IAU 2006/2000A, CIO based, using X, Y series" in the SOFA Tools for Earth Attitude [1].

    Args:
        time:  A lib.time Time-object

    Returns:
        numpy.array: An array of 3x3 transformation matrices
    """
    return R3(-E(time)) @ R2(-d(time)) @ R3(E(time)) @ R3(s(time))


@lru_cache()
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

    return R3(-ERA(time))


@lru_cache()
def W(time):
    """Transformation matrix for the polar motion

    The transformation matrix is described in the IERS Conventions [2], section 5.4.1. The implementation is based on
    section 5.6 "IAU 2006/2000A, CIO based, using X, Y series" in the SOFA Tools for Earth Attitude [1].

    Args:
        time:  A lib.time Time-object

    Returns:
        numpy.array: An array of 3x3 transformation matrices
    """
    return R3(-s_prime(time)) @ R2(xp(time)) @ R1(yp(time))


#
# Rotation angles
#
@lru_cache()
def s(time):
    """CIO locator
    """
    return sofa.vectorized_s06(time)


@lru_cache()
def s_prime(time):
    """TIO locator
    """
    return sofa.vectorized_sp00(time)


@lru_cache()
def ERA(time):
    """Earth rotation angle
    """
    return sofa.vectorized_era00(time)


@lru_cache()
def d(time):
    """Polar angle of CIP from positive z-axis

    Spherical coordinate with radius = 1
    """
    return np.arccos(Z(time))


@lru_cache()
def E(time):
    """Azimuth angle of CIP in xy-plane

    Spherical coordinate with radius = 1
    """
    return np.arctan2(Y(time), X(time))


@lru_cache()
def xp(time):
    """X coordinate of the CIP in ITRS
    """
    eop = apriori.get("eop", time=time)
    return eop.x * Unit.arcsec2rad


@lru_cache()
def yp(time):
    """Y coordinate of the CIP in ITRS
    """
    eop = apriori.get("eop", time=time)
    return eop.y * Unit.arcsec2rad


@lru_cache()
def X(time):
    """Total X coordinate of CIP

    Returns:
        model + celestial pole offset
    """
    eop = apriori.get("eop", time=time)
    return sofa.X_model(time) + eop.dx * Unit.arcsec2rad


@lru_cache()
def Y(time):
    """Total Y coordinate of CIP

    Returns:
        model + celestial pole offset
    """
    eop = apriori.get("eop", time=time)
    return sofa.Y_model(time) + eop.dy * Unit.arcsec2rad


@lru_cache()
def Z(time):
    """Total Z coordinate of CIP

    Derived from X and Y with assumption of unit sphere with

    Returns:
        model + celestial pole offset
    """
    return np.sqrt(1 - (X(time) ** 2 + (Y(time) ** 2)))


#
# Rotation matricies differentiated with regards to earth orientation parameters
#


@lru_cache()
def dW_dxp(time):
    """Derivative of transformation matrix for the polar motion with regards to the CIP (Celestial Intermediate Pole)
    in TRF along the Greenwich meridian x_p.

    This is done according to equations (2.31) in Teke :cite:`teke2011` which is the analytical partial derivative of
    equation (5.3) in :cite:'iers2010'.
    """
    return R3(-s_prime(time)) @ dR2(xp(time)) @ R1(yp(time))


@lru_cache()
def dW_dyp(time):
    """Derivative of transformation matrix for the polar motion with regards to the CIP in ITRS.

    Analytical partial derivative of equation (5.3) in :cite:'iers2010'.
    """
    return R3(-s_prime(time)) @ R2(xp(time)) @ dR1(yp(time))


@lru_cache()
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

    return -dR3(-ERA(time)) * constant.omega


@lru_cache()
def dQ_dX(time):
    """Derivative of transformation matrix for nutation/presession with regards to the X coordinate of CIP in GCRS
    """
    # Rotation matrices
    R3_E = R3(E(time))
    R3_s = R3(s(time))
    R2_md = R2(-d(time))
    R3_mE = R3(-E(time))
    dR3_s = dR3(s(time))
    dR3_E = dR3(E(time))
    dR3_mE = dR3(-E(time))
    dR2_md = dR2(-d(time))

    return (
        dR3_mE @ R2_md @ R3_E @ R3_s * (-dE_dX(time))
        + R3_mE @ dR2_md @ R3_E @ R3_s * (-dd_dX(time))
        + R3_mE @ R2_md @ dR3_E @ R3_s * (dE_dX(time))
        + R3_mE @ R2_md @ R3_E @ dR3_s * (ds_dX(time))
    )


@lru_cache()
def dQ_dY(time):
    """Derivative of transformation matrix for nutation/presession with regards to the Y coordinate of CIP in GCRS
    """
    # Rotation matrices
    R3_E = R3(E(time))
    R3_s = R3(s(time))
    R2_md = R2(-d(time))
    R3_mE = R3(-E(time))
    dR3_s = dR3(s(time))
    dR3_E = dR3(E(time))
    dR3_mE = dR3(-E(time))
    dR2_md = dR2(-d(time))

    return (
        dR3_mE @ R2_md @ R3_E @ R3_s * (-dE_dY(time))
        + R3_mE @ dR2_md @ R3_E @ R3_s * (-dd_dY(time))
        + R3_mE @ R2_md @ dR3_E @ R3_s * (dE_dY(time))
        + R3_mE @ R2_md @ R3_E @ dR3_s * (ds_dY(time))
    )


#
# Rotation angles differentiated with regards to earth orientation parameters
#


@lru_cache()
def dE_dX(time):
    """Derivative of azimuth angle of CIP in GCRS with regards to the X coordinate of CIP in GCRS
    """
    return (-Y(time) / (X(time) ** 2 + Y(time) ** 2))[:, None, None]


@lru_cache()
def ds_dX(time):
    """Derivative of CIO locator with regards the X coordinate of CIP in GCRS
    """
    return (-Y(time) / 2)[:, None, None]


@lru_cache()
def dd_dX(time):
    """Derivative of polar angle of CIP in GCRS with regards to the X coordinate of CIP in GCRS
    """
    return (X(time) / (Z(time) * np.sqrt(X(time) ** 2 + Y(time) ** 2)))[:, None, None]


@lru_cache()
def dE_dY(time):
    """Derivative of azimuth angle of CIP in GCRS with regards to the Y coordinate of CIP in GCRS
    """
    return (X(time) / (X(time) ** 2 + Y(time) ** 2))[:, None, None]


@lru_cache()
def ds_dY(time):
    """Derivative of CIO locator with regards the X coordinate of CIP in GCRS
    """
    return (-X(time) / 2)[:, None, None]


@lru_cache()
def dd_dY(time):
    """Derivative of polar angle of CIP in GCRS with regards to the Y coordinate of CIP in GCRS
    """
    return (Y(time) / (Z(time) * np.sqrt(X(time) ** 2 + Y(time) ** 2)))[:, None, None]
