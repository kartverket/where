# cython: profile=True
"""Calculates the force acting on the satellite from the drag from particles in the Earth atmosphere

Description:

This model calculates the drag force from particles in the Earth atmosphere on the satellite.

References:

   [1] O. Montenbruck and E. Gill: Satellite Orbits, Springer, 2000.
"""

# Standard library imports
import numpy as np
import math

# Where imports
from where import apriori

cdef double drag_coefficient
cdef double area
cdef double mass
cdef int drag_idx
cdef double[:, :, :] g2i


def register_entry_point():
    """Register entry points for setup and later calls."""
    return dict(setup=drag_setup, call=drag)


def drag_setup(
        rundate, force_parameters, sat_name, time_grid, epochs, body_pos_gcrs, body_pos_itrs, bodies, gcrs2itrs
):
    """Set up module variables used later during calculation.

    Args:
        rundate:           Time of integration start
        force_parameters:  Dict of parameters to be estimated
        sat_name:          Name of satellite
        time_grid:         Table of times in seconds since rundate, in utc.
        epochs:            time_grid converted to Time objects, in utc.
        body_pos_gcrs:     The positions of the bodies in the solar system in GCRS.
        body_pos_itrs:     The positions of the bodies in the solar system in ITRS.
        bodies:            List of bodies
        gcrs2itrs:         List of transformation matrices, one for each time in epochs.
    """
    global drag_coefficient, area, mass, drag_idx, g2i

    sat = apriori.get_satellite(sat_name)
    drag_idx = -1
    drag_coefficient = 0.0
    if "drag_coefficient" in force_parameters:
        drag_idx = list(force_parameters.keys()).index("drag_coefficient")
        drag_coefficient = force_parameters["drag_coefficient"]
    else:
        drag_coefficient = sat.drag_coefficient

    area = sat.area
    mass = sat.mass

    g2i = gcrs2itrs


def drag(double[:] sat_vel_itrs, int num_param, int current_step, **_not_used):
    """Compute drag force on satellite

    Args:
        sat_vel_itrs:      Satellite velocity in ITRS.
        num_param:         Number of parameters to be estimated.
        current_step:      Int, step number of current step of integrator.
        _not_used:         Unused variables.

    Returns:
        Acceleration and equation for state transition matrix due
        to drag from particles in the atmosphere of the Earth.
    """
    cdef double sat_speed_itrs, sat_speed_squared
    cdef double v1, v2, v3
    cdef double[:] acc, dacc_dp, acc_gcrs, dacc_dp_gcrs
    cdef double[:, :] trans_itrs, trans_gcrs, trans, sens
    cdef int i

    acc = np.zeros(3)
    acc_gcrs = np.zeros(3)
    dacc_dp = np.zeros(3)
    v1, v2, v3 = sat_vel_itrs
    sat_speed_squared = v1**2 + v2**2 + v3**2
    sat_speed_itrs = math.sqrt(sat_speed_squared)
    sens = np.zeros((3, num_param))

    # Equation 3.97 from Montenbruck and Gill [2], with drag_coefficient being the product of the
    # actual drag_coefficient and the atmospheric density, both unknown. Estimate their product.
    acc[0] = -0.5 * drag_coefficient * area * sat_speed_itrs * v1 / mass
    acc[1] = -0.5 * drag_coefficient * area * sat_speed_itrs * v2 / mass
    acc[2] = -0.5 * drag_coefficient * area * sat_speed_itrs * v3 / mass

    # The derivative of the acceleration with respect to the relevant parameter
    dacc_dp[0] = -0.5 * area * sat_speed_itrs * v1 / mass
    dacc_dp[1] = -0.5 * area * sat_speed_itrs * v2 / mass
    dacc_dp[2] = -0.5 * area * sat_speed_itrs * v3 / mass

    # Equation for state transition matrix
    trans_itrs = -np.array(((2 * v1**2 + v2**2 + v3**2, v1*v2, v1*v3),
                            (v1*v2, 2 * v2**2 + v1**2 + v3**2, v2*v3),
                            (v1*v3, v2*v3, 2 * v3**2 + v1**2 + v2**2)
                           )) * ((0.5 * drag_coefficient * area) / (mass * sat_speed_itrs))
    # Transform to space fixed system before returning
    gcrs2itrs = g2i[current_step]
    acc_gcrs = np.dot(gcrs2itrs.T, acc)
    trans_gcrs = np.dot(gcrs2itrs.T, np.dot(trans_itrs, gcrs2itrs))
    dacc_dp_gcrs = np.dot(gcrs2itrs.T, dacc_dp)

    trans = np.hstack((np.zeros((3, 3)), trans_gcrs))

    if not drag_idx == -1:
        for j in range(0, 3):
            sens[j, drag_idx] = dacc_dp_gcrs[j]
    return (acc_gcrs, trans, sens)
