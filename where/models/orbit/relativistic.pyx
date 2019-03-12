# cython: profile=True
"""Calculates the force acting on the satellite from relativistic effects

Description:

Calculates the force acting on the satellite from relativistic effects, following section 3.7.3 in [1].

References:

   [1] O. Montenbruck and E. Gill: Satellite Orbits, Springer, 2000.

"""

# Standard library imports
import sys

# External library imports
import numpy as np

# Where imports
from midgard.math.constant import constant
from where.lib import log

cdef double GM, c
cdef int num_param
cdef double [:, :, :] g2i


def register_entry_point():
    """Register entry points for setup and later calls."""
    return dict(setup=relativistic_setup, call=relativistic)


def relativistic_setup(
        rundate, force_parameters, sat_name, time_grid, epochs, body_pos_gcrs, body_pos_itrs, bodies, gcrs2itrs
):
    """Set up module variables used later during calculation.

    Args:
        rundate:           Time of integration start.
        force_parameters:  Dict of parameters to be estimated.
        sat_name:          Name of satellite.
        time_grid:         Table of times in seconds since rundate, in utc.
        epochs:            time_grid converted to Time objects, in utc.
        body_pos_gcrs:     The positions of the bodies in the solar system in GCRS.
        body_pos_itrs:     The positions of the bodies in the solar system in ITRS.
        bodies:            List of bodies.
        gcrs2itrs:         List of transformation matrices, one for each time in epochs.
    """
    global GM, c, num_param
    global g2i

    # Set gravitational constant and speed of light
    GM = constant.get("GM", source="egm_2008")
    c = constant.get("c")

    # Number of parameters to be estimated
    num_param = len(force_parameters)

    if not (sat_name == "lageos1" or sat_name == "lageos2"):
        log.fatal(f"Unknown relativistic effects for satellite {sat_name}")
        sys.exit(0)
    g2i = gcrs2itrs


def relativistic(sat_pos_itrs, sat_vel_itrs, int current_step, **_not_used):
    """Compute relativistic gravitational force on satellite

    Args:
        sat_pos_itrs:     Satellite position in ITRS.
        sat_vel_itrs:     Satellite velocity in ITRS.
        current_step:     Int, step number of current step of integrator.
        _not_used:        Unused variables.

    Returns:
        Acceleration and equation for state transition matrix due to relativistic effects.
    """
    cdef int i
    cdef double[:] acc = np.zeros(3)

    # Equation 3.147 from Montenbruck and Gill [1].
    # Assume circular orbit.

    acc = -GM * sat_pos_itrs / np.linalg.norm(sat_pos_itrs)**3 * 3 * sat_vel_itrs**2 / c**2
    gcrs2itrs = g2i[current_step]
    # Transform to space fixed system before returning
    acc = np.dot(gcrs2itrs.T, acc)

    # Assume negligible state transition matrix
    trans = np.zeros((3, 6))

    # No parameters to be estimated
    sens = np.zeros((3, num_param))

    return (acc, trans, sens)
