"""Library for basic rotation matrices


Description:
------------
Creates rotation matrices for rotation around the axes of a right handed Cartesian coordinate system and
their derivatives.

For instance, for an XYZ-system, R1 returns a rotation matrix around the x-axis and for an ENU-system, R1 returns a
rotation matrix around the east-axis. dR1 returns the derivative of the R1 matrix with respect to the rotation
angle. All functions are vectorized, so that one rotation matrix is returned per input angle.

Example:

>>> from where.lib import rotation
>>> rotation.R1([0, 1])
array([[[ 1.        ,  0.        ,  0.        ],
        [ 0.        ,  1.        ,  0.        ],
        [ 0.        , -0.        ,  1.        ]],

       [[ 1.        ,  0.        ,  0.        ],
        [ 0.        ,  0.54030231,  0.84147098],
        [ 0.        , -0.84147098,  0.54030231]]])




$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""

import numpy as np


def R1(angle):
    """Rotation matrix around the first axis

    Args:
        angle:  Scalar, list or numpy array of angles in radians.

    Returns:
        Numpy array:   Rotation matrix or array of rotation matrices.
    """
    zero, one = _zero(angle), _one(angle)
    cosA, sinA = np.cos(angle), np.sin(angle)
    return _roll_axes(np.array([[one, zero, zero], [zero, cosA, sinA], [zero, -sinA, cosA]]))


def R2(angle):
    """Rotation matrix around the second axis

    Args:
        angle:  Scalar, list or numpy array of angles in radians.

    Returns:
        Numpy array:   Rotation matrix or array of rotation matrices.
    """
    zero, one = _zero(angle), _one(angle)
    cosA, sinA = np.cos(angle), np.sin(angle)
    return _roll_axes(np.array([[cosA, zero, -sinA], [zero, one, zero], [sinA, zero, cosA]]))


def R3(angle):
    """Rotation matrix around the third axis

    Args:
        angle:  Scalar, list or numpy array of angles in radians.

    Returns:
        Numpy array:   Rotation matrix or array of rotation matrices.
    """
    zero, one = _zero(angle), _one(angle)
    cosA, sinA = np.cos(angle), np.sin(angle)
    return _roll_axes(np.array([[cosA, sinA, zero], [-sinA, cosA, zero], [zero, zero, one]]))


def dR1(angle):
    """Derivative of a rotation matrix around the first axis with respect to the rotation angle.

    Args:
        angle:  Scalar, list or numpy array of angles in radians.

    Returns:
        Numpy array:   Rotation matrix or array of rotation matrices.
    """
    zero = _zero(angle)
    cosA, sinA = np.cos(angle), np.sin(angle)
    return _roll_axes(np.array([[zero, zero, zero], [zero, -sinA, cosA], [zero, -cosA, -sinA]]))


def dR2(angle):
    """Derivative of a rotation matrix around the second axis with respect to the rotation angle

    Args:
        angle:  Scalar, list or numpy array of angles in radians.

    Returns:
        Numpy array:   Rotation matrix or array of rotation matrices.
    """
    zero = _zero(angle)
    cosA, sinA = np.cos(angle), np.sin(angle)
    return _roll_axes(np.array([[-sinA, zero, -cosA], [zero, zero, zero], [cosA, zero, -sinA]]))


def dR3(angle):
    """Derivative of a rotation matrix around the third axis with respect to the rotation angle

    Args:
        angle:  Scalar, list or numpy array of angles in radians.

    Returns:
        Numpy array:   Rotation matrix or array of rotation matrices.
    """
    zero = _zero(angle)
    cosA, sinA = np.cos(angle), np.sin(angle)
    return _roll_axes(np.array([[-sinA, cosA, zero], [-cosA, -sinA, zero], [zero, zero, zero]]))


def enu2trf(lat, lon):
    """Rotation matrix for rotating an ENU coordinate system to an earth oriented one

    See for instance http://www.navipedia.net/index.php/Transformations_between_ECEF_and_ENU_coordinates
    This is equal to doing::

        R3(-(np.pi/2 + lon)) @ R1(-(np.pi/2 - lat))

    Args:
        lat (Float or Array):   Latitude of origin of ENU coordinate system.
        lon (Float or Array):   Longitude of origin of ENU coordinate system.

    Returns:
        Numpy array:   Rotation matrix or array of rotation matrices.
    """
    zero = _zero(lat)
    coslat, coslon, sinlat, sinlon = np.cos(lat), np.cos(lon), np.sin(lat), np.sin(lon)
    return _roll_axes(
        np.array(
            [
                [-sinlon, -coslon * sinlat, coslon * coslat],
                [coslon, -sinlon * sinlat, sinlon * coslat],
                [zero, coslat, sinlat],
            ]
        )
    )


def trf2enu(lat, lon):
    """Rotation matrix for rotating an earth oriented coordinate system to an ENU one

    See for instance http://www.navipedia.net/index.php/Transformations_between_ECEF_and_ENU_coordinates
    This is equal to doing::

        R1(np.pi/2 - lat) @ R3(np.pi/2 + lon)

    Args:
        lat (Float or Array):   Latitude of origin of ENU coordinate system.
        lon (Float or Array):   Longitude of origin of ENU coordinate system.

    Returns:
        Numpy array:   Rotation matrix or array of rotation matrices.
    """
    zero = _zero(lat)
    coslat, coslon, sinlat, sinlon = np.cos(lat), np.cos(lon), np.sin(lat), np.sin(lon)
    return _roll_axes(
        np.array(
            [
                [-sinlon, coslon, zero],
                [-sinlat * coslon, -sinlat * sinlon, coslat],
                [coslat * coslon, coslat * sinlon, sinlat],
            ]
        )
    )


def _roll_axes(mat):
    """Move the axes of an array of 2-d matrices properly

    Roll the first two dimensions to the end, so that indexing works as expected.

    Args:
        mat (numpy array):  Array of 2-d matrices, can be a numpy array of any dimension.

    Returns:
        Numpy array:  The same array as mat, but with the first two dimensions rolled to the end.
    """
    return mat if mat.ndim < 3 else mat.transpose(np.roll(np.arange(mat.ndim), -2))


def _zero(angle):
    """Returns a scalar or array of zeros with the same shape as angle

    Args:
        angle:   Scalar, list or numpy array.

    Returns:
        Scalar or numpy array:  Zero-scalar or array.
    """
    try:
        return np.zeros(angle.shape)  # angle is numpy array
    except AttributeError:
        try:
            return np.zeros(len(angle))  # angle is list or other iterable
        except TypeError:
            return 0  # angle is scalar


def _one(angle):
    """Returns a scalar or array of ones with the same shape as angle

    Args:
        angle:   Scalar, list or numpy array.

    Returns:
        Scalar or numpy array:  One-scalar or array.
    """
    try:
        return np.ones(angle.shape)  # angle is numpy array
    except AttributeError:
        try:
            return np.ones(len(angle))  # angle is list or other iterable
        except TypeError:
            return 1  # angle is scalar
