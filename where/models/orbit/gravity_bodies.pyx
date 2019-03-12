# cython: profile=True
"""Calculates the force on the satellite from the gravity field of the Sun, Moon and planets

Description:

Calculates the gravitational force from the bodies, assuming they are point masses, following
Montenbruck and Gill [1].

References:
[1] Montenbruck, Oliver and Gill, Eberhard, Satellite Orbits,
    Springer Verlag, 2000.
"""

# External library imports
import numpy as np
import cython
import math

# Where imports
from where import apriori
from midgard.math.constant import constant
from where.lib import config

# Name of model
MODEL = __name__.split(".")[-1]

# Module variables set during setup
cdef double[:] GM_bodies
cdef double[:, :, :] body_pos
cdef int num_bodies


def register_entry_point():
    """Register entry points for setup and later calls."""
    return dict(setup=gravity_bodies_setup, call=gravity_bodies)


def gravity_bodies_setup(
        rundate, force_parameters, sat_name, time_grid, epochs, body_pos_gcrs, body_pos_itrs, bodies, gcrs2itrs
):
    """Set up module variables used later during calculation.

    Args:
        rundate:            Time of integration start.
        force_parameters:   Dict of parameters to be estimated.
        sat_name:           Name of satellite.
        time_grid:          Table of times in seconds since rundate, in utc.
        epochs:             time_grid converted to Time objects, in utc.
        body_pos_gcrs:      The positions of the bodies in the solar system in GCRS.
        body_pos_itrs:      The positions of the bodies in the solar system in ITRS.
        bodies:             List of bodies.
        gcrs2itrs:         List of transformation matrices, one for each time in epochs.
    """
    global num_bodies, GM_bodies, body_pos
    cdef int idx, i
    body_pos = body_pos_gcrs
    num_bodies = len(bodies)
    GM_bodies = np.zeros(num_bodies)

    for idx in range(num_bodies):
        GM_bodies[idx] = constant.get("GM_" + bodies[idx])


@cython.boundscheck(False)
@cython.wraparound(False)
def gravity_bodies(double[:] sat_pos_gcrs, force_parameters, int current_step, **_not_used):
    """Compute force on satellite from the gravity field of the bodies

    The force is calculated assuming that the Moon is a point mass. For the acceleration we have from equation (3.37):

    \f[ \ddot{\vec r} = GM \cdot \left( \frac{\vec s - \vec r}{| \vec s -
    \vec r |^3} - \frac{\vec s}{| \vec s |^3} \right) . \f]

    Here \f$ \vec r \f$ is the geocentric position vector of the satellite, while \f$ \vec s \f$ is the geocentric
    position vector of the Moon.  For the transition matrix we use equation (7.75):

    \f[ \frac{\partial \ddot{\vec r}}{\partial \vec r} = -GM \cdot \left(
    \frac{1}{| \vec r - \vec s |^3} I_{3 \times 3} - 3 (\vec r - \vec s)
    \frac{(\vec r - \vec s)^T}{| \vec r - \vec s |^5} \right) . \f]

    Args:
        sat_pos_gcrs:      Satellite position as 3-vector in GCRS.
        force_parameters:  ordered dict of parameters to be estimated
        _not_used:         Unused variables.

    Returns:
        Acceleration and transition matrix due to gravity from the bodies in GCRS.
    """
    cdef double[:] acc, sat_body_norm
    cdef double[:, :] trans, sat_body_vec
    cdef int idx, i

    acc = np.zeros(3)
    trans = np.zeros((3, 3))
    sat_body_vec = np.zeros((num_bodies, 3))
    sat_body_norm = np.zeros(num_bodies)
    body_pos_norm = np.zeros(num_bodies)
    body_pos_temp = np.zeros((num_bodies, 3))

    for idx in range(num_bodies):
        body_pos_temp[idx, :] = body_pos[idx, current_step]
    for idx in range(num_bodies):
        for i in range(0, 3):
            sat_body_vec[idx, i] = body_pos_temp[idx, i] - sat_pos_gcrs[i]
            sat_body_norm[idx] += sat_body_vec[idx, i]**2
            body_pos_norm[idx] += body_pos_temp[idx, i]**2

        sat_body_norm[idx] = math.sqrt(sat_body_norm[idx])
        body_pos_norm[idx] = math.sqrt(body_pos_norm[idx])

        for i in range(0, 3):
            acc[i] += GM_bodies[idx] * (sat_body_vec[idx, i] / sat_body_norm[idx]**3
                                        - body_pos_temp[idx, i] / body_pos_norm[idx]**3)
            for j in range(0, 3):
                trans[i, j] += -GM_bodies[idx] * ((i == j) / sat_body_norm[idx]**3 - 3 * sat_body_vec[idx, i] *
                                                  sat_body_vec[idx, j] / sat_body_norm[idx]**5)

    return acc, np.hstack((trans, np.zeros((3, 3)))), np.zeros((3, len(force_parameters)))
