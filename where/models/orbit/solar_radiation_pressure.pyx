# cython: profile=True
"""Calculates the satellite position at observation epochs

Description:

Calculates the acceleration on the satellite caused by the solar radiation pressure. The calculation follows
Montenbruck and Gill [1].

References:
[1] Montenbruck, Oliver and Gill, Eberhard, Satellite Orbits, Springer Verlag, 2000.
"""

# Standard library imports
cimport libc.math
import numpy as np

# Where imports
from where import apriori
from midgard.math.constant import constant
from where.lib import plugins

cdef double radiation_pressure_coefficient
cdef double c
cdef double[:, :] sun_pos
cdef int rad_idx
cdef satellite
cdef sun_flux
cdef double[:, :] ident


def register_entry_point():
    """Register entry points for setup and later calls."""
    return dict(call=solar_radiation_pressure, setup=solar_setup)


def solar_setup(
        rundate, force_parameters, sat_name, time_grid, epochs, body_pos_gcrs, body_pos_itrs, bodies, gcrs2itrs
):
    """
    Args:
        rundate:           Rundate, used as a reference date.
        force_parameters:  Dict of parameters to be estimated.
        sat_name:          Name of satellite.
        time_grid:         Table of times in seconds since rundate, in utc.
        epochs:            time_grid converted to Time objects, in utc.
        body_pos_gcrs:     The positions of the bodies in the solar system in GCRS.
        body_pos_itrs:     The positions of the bodies in the solar system in ITRS.
        bodies:            List of bodies.
        gcrs2itrs:         List of transformation matrices, one for each time in epochs.
    """
    global satellite, rad_idx, radiation_pressure_coefficient
    global c
    global sun_flux
    global sun_pos
    global ident

    ident = np.eye(3)
    satellite = apriori.get_satellite(sat_name)

    # Old code using data for the solar flux 
    # Not in use at the moment
    # flux_table = apriori.get("solar_flux", rundate=rundate)
    # sun_flux = flux_table(time_grid)

    rad_idx = -1
    if "radiation_pressure_coefficient" in force_parameters:
        rad_idx = list(force_parameters.keys()).index("radiation_pressure_coefficient")
        radiation_pressure_coefficient = force_parameters["radiation_pressure_coefficient"]
    else:
        radiation_pressure_coefficient = satellite.radiation_pressure_coefficient

    c = constant.get("c")
    flux = constant.get("S", source="book")
    
    sun_flux = np.repeat(flux, len(time_grid))
    idx = bodies.index("sun")
    sun_pos = body_pos_gcrs[idx, :, :]


def solar_radiation_pressure(double[:] sat_pos_gcrs, str sat_name, int num_param, int current_step, **_not_used):
    """Compute force on satellite from the solar radiation pressure

    The force on a satellite caused by the solar radiation pressure is
    described in chapter 3.4 of Montenbruck and Gill [1].

    Args:
        sat_pos_gcrs:                    Satellite position as a 3-vector in GCRS.
        num_param:                       Number of parameters to be estimated.
        sat_name:                        Name of satellite.
        _not_used:                       Unused variables.

    Returns:
        Acceleration and equation for state transition matrix due
        to solar radiation pressure acting on the satellite.

    """
    global ident
    cdef double[:] sat_sun_vec = np.zeros(3)
    cdef int i
    cdef double[:, :] trans = np.zeros((3, 6))

    for i in range(0, 3):
        sat_sun_vec[i] = sun_pos[current_step, i] - sat_pos_gcrs[i]
    # Calculate shadow function based on apparent radius of sun and earth
    cdef double sat_sun_norm = np.linalg.norm(sat_sun_vec)
    cdef double sat_earth_norm = np.linalg.norm(sat_pos_gcrs)
    cdef double sun_rad = libc.math.asin(constant.R_sun / sat_sun_norm)
    cdef double earth_rad = libc.math.asin(constant.a / sat_earth_norm)
    cdef double sun_earth_sep = libc.math.acos(-np.dot(sat_pos_gcrs, sat_sun_vec) / (sat_earth_norm * sat_sun_norm))
    
    # Acceleration of satellite due to solar radiation pressure
    # Equation (3.75) in Montenbruck [1]
    cdef double[:] acc = np.zeros(3)
    for i in range(0, 3):
        acc[i] = (
            -radiation_pressure_coefficient * (sun_flux[current_step] / c) * (satellite.area / satellite.mass)
            * sat_sun_vec[i] / sat_sun_norm
        )

    # Equation for state transition matrix
    # Equation (7.77) in Montenbruck [1]
    # Removed the AU*2 / sat_sun_norm**2 factor, since we use observed flux values instead of average values
    for i in range(0, 3):
        for j in range(0, 3):
            trans[i, j] = (
                radiation_pressure_coefficient * (sun_flux[current_step] / c) * (satellite.area / satellite.mass)
                / sat_sun_norm * ((i == j) - 3 * sat_sun_vec[i] * sat_sun_vec[j] / sat_sun_norm**2)
            )

    # Derivative of acceleration with respect to radiation pressure coefficient
    cdef double[:] dacc_dp = np.zeros(3)
    for i in range(0, 3):
        dacc_dp[i] = -(sun_flux[current_step] / c) * (satellite.area / satellite.mass) * sat_sun_vec[i] / sat_sun_norm

    cdef double[:, :] sens = np.zeros((3, num_param))

    if not rad_idx == -1:
        for i in range(0, 3):
            sens[i, rad_idx] = dacc_dp[i]
    # The easy cases: Satellite in total shadow or full sunlight
    if sun_earth_sep < earth_rad - sun_rad:
        return np.zeros(3), np.zeros((3, 6)), np.zeros((3, num_param))
    elif earth_rad + sun_rad <= sun_earth_sep:
        return acc, trans, sens

    # Satellite in partial shadow
    x = (sun_earth_sep**2 + sun_rad**2 - earth_rad**2) / (2 * sun_earth_sep)
    y = libc.math.sqrt(sun_rad**2 - x**2)
    shadow_area = (
        sun_rad**2 * libc.math.acos(x / sun_rad) + earth_rad**2 * libc.math.acos((sun_earth_sep - x) / earth_rad)
        - sun_earth_sep * y
    )
    scaling_factor = 1 - shadow_area / (libc.math.pi * sun_rad**2)

    for i in range(0, 3):
        acc[i] *= scaling_factor
        for j in range(0, 6):
            trans[i, j] *= scaling_factor
        for k in range(0, num_param):
            sens[i, k] *= scaling_factor
    return (acc, trans, sens)
