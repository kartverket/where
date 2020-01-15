# cython: profile=True
"""Calculates the force acting on the satellite from the infrared emissivity and optical albedo of the Earth

Description:

This model calculates the indirect solar radiation pressure (both infrared
and optical) from the Sun reflected by Earth.

References:

   [1] Earth Radiation Pressure Effects on Satellites, P.C. Knocke, J.C. Ries, B. D. Tapley,
       Proceedings of the AIAA/AAS Astrodynamics Conference, pp 577-587 (1988).
"""
# Standard library imports
import math

# External library imports
import numpy as np

# Where imports
from where import apriori
from midgard.math.constant import constant
from where.ext import sofa
from where.data.time import Time, TimeDelta

cdef double c, earth_radius, AU
cdef double area, mass
cdef double radiation_pressure_coefficient
cdef double[:] omega
cdef double[:] solar_flux_over_c
cdef double[:, :] sun_pos_unit
cdef double[:, :, :] g2i

def register_entry_point():
    """Register entry points for setup and later calls."""
    return dict(setup=knocke_setup, call=indirect_radiation)


def knocke_setup(
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
    global c, earth_radius, AU
    global area, mass
    global solar_flux_over_c
    global radiation_pressure_coefficient
    global omega
    global sun_pos_unit
    global g2i

    cdef int i, j
    cdef double days = 0
    omega = np.zeros(len(epochs))
    solar_flux_over_c = np.zeros(len(epochs))
    sun_pos_unit = np.zeros((len(epochs), 3))

    sat = apriori.get_satellite(sat_name)
    earth_radius = constant.get("a", source="egm_2008")
    AU = constant.get("AU", source="web")
    c = constant.get("c")
    flux = constant.get("S")
    area = sat.area
    mass = sat.mass

    # Code for using flux data, not in use. 
    # Using only constant value at the moment.
    # flux_table = apriori.get("solar_flux", rundate=rundate)
    # sun_flux = flux_table(time_grid)

    sun_flux = np.repeat(flux, len(epochs))
    if "radiation_pressure_coefficient" in force_parameters:
        radiation_pressure_coefficient = force_parameters["radiation_pressure_coefficient"]
    else:
        radiation_pressure_coefficient = sat.radiation_pressure_coefficient

    sun_pos_itrs = body_pos_itrs[bodies.index("sun"), :, :]

    earth_sun_distance_in_au = np.linalg.norm(sun_pos_itrs, axis=1) / AU
    sun_pos_norm = np.linalg.norm(sun_pos_itrs, axis=1)

    for i in range(0, len(epochs)):
        # Days since reference epoch Dec 22 1981
        days = (epochs.tt[i] - Time("1981-12-22", scale="tt", fmt="date")).jd
        omega[i] = 2 * math.pi * days / 365.25

        # This is E_S/c in equation 27 in Knocke:
        solar_flux_over_c[i] = sun_flux[i] / (c * earth_sun_distance_in_au[i]**2)
        for j in range(0, 3):
            sun_pos_unit[i, j] = sun_pos_itrs[i, j] / sun_pos_norm[i]
    g2i = gcrs2itrs


def indirect_radiation(double[:] sat_pos_itrs, int num_param, int current_step, **_not_used):
    """Compute the force acting on the satellite from the radiation pressure of the Earth, following Knocke [1]

    Args:
        sat_pos_itrs:             Satellite position in ITRS.
        num_param:                int, number of parameters to be estimated.
        current_step:             Current step number of integrator.
        _not_used:                Unused variables.

    Returns:
        Acceleration and equation for state transition matrix due to infrared
        emissivity and reflection of visible light of the Earth in GCRS.
    """
    acc = np.zeros(3)
    trans = np.zeros((3, 3))
    cdef int i, j
    #cdef double latitude, longitude
    #cdef double sin_latitude, cos_latitude
    #cdef double ka, ke

    latitude = 0
    longitude = 0
    sin_latitude = 0
    cos_latitude = 0
    omega_now = omega[current_step]
    cos_omega = math.cos(omega_now)

    # Create a grid of n deg x n deg surface elements covering the earth.
    # n=steplength:
    steplength = 15

    for latitude in range(-90, 90, steplength):
        # Consider mid-point of surface element:
        latitude += 0.5 * steplength
        latitude = math.radians(latitude)
        sin_latitude = math.sin(latitude)
        cos_latitude = math.cos(latitude)
        ka = knocke_a(sin_latitude, cos_omega)
        ke = knocke_e(sin_latitude, cos_omega)
        for longitude in range(0, 360, steplength):
            # Consider mid-point of surface element:
            longitude += 0.5 * steplength
            longitude = math.radians(longitude)

            surface_elem_vector, j = sofa.iau_gd2gc(2, longitude, latitude, 0)
            surface_elem_sat_vector = sat_pos_itrs - surface_elem_vector
            # To be divided by norms later
            cos_sat_zenith_angle = np.dot(surface_elem_sat_vector, surface_elem_vector)

            if cos_sat_zenith_angle < 0:
                #Satellite not visible from surface element
                continue

            surface_elem_distance = np.linalg.norm(surface_elem_vector)
            surface_elem_sat_distance = np.linalg.norm(surface_elem_sat_vector)
            cos_sat_zenith_angle *= 1 / (surface_elem_sat_distance * surface_elem_distance)

            # width of surface element
            width = 2 * earth_radius * steplength * math.pi * cos_latitude / 360

            # Height of surface element
            height = 2 * steplength * math.pi * earth_radius / 360
            area_surface_elem = height * width
            effective_area = area_surface_elem * cos_sat_zenith_angle


            # Cos of zenith angle of the Sun.
            cos_theta = np.dot(surface_elem_vector / surface_elem_distance, sun_pos_unit[current_step])
            if cos_theta < 0:
            # Night-time, no reflection of sunlight
            # Shadow function is zero.
                cos_theta = 0

            # Satellite acceleration due to earth emissivity:
            acc += ((ka * cos_theta + 0.25 * ke) * effective_area * surface_elem_sat_vector
                    / surface_elem_sat_distance**3
                   )
            # Equation for state transition matrix
            v = ((ka * cos_theta + 0.25 * ke) * area_surface_elem * surface_elem_sat_distance**(-6)
                 * surface_elem_distance**(-1)
                )

            r1, r2, r3 = surface_elem_sat_vector[0: 3]
            rs1, rs2, rs3 = surface_elem_vector[0: 3]
            trans += v * np.array(
                      [[2 * rs1 * r1 * (r2**2 + r3**2 - r1**2), -8 * r1 * 2 * rs1 * r2, -8 * r1**2 * rs1 * r3],
                       [-8 * r1 * r2**2 * rs2, 2 * rs2 * r2 * (r1**2 + r3**2 - r2**2), -8 * r2**2 * rs2 * r3],
                       [-8 * r1 * r3**2 * rs3, -8 * r2 * r3**2 * rs3, 2 * rs3 * r3 * (r1**2 + r2**2 - r3**2)]
                      ])
    for i in range(0, 3):
        acc[i] *= radiation_pressure_coefficient * solar_flux_over_c[current_step] * (area / mass) * (1 / math.pi)
        for j in range(0, 3):
            trans[i, j] *= (
                radiation_pressure_coefficient * solar_flux_over_c[current_step] * (area / mass) * (1 / math.pi)
            )

    # Transform to space fixed system before returning
    gcrs2itrs = g2i[current_step]
    acc = np.dot(gcrs2itrs.T, acc)
    trans = np.dot(gcrs2itrs.T, np.dot(trans, gcrs2itrs))
    #dacc_dp = np.dot(gcrs2itrs.T, dacc_dp)

    trans = np.hstack((trans, np.zeros((3, 3))))
    sens = np.zeros((3, num_param))

    return (acc, trans, sens)


def knocke_e(sin_latitude, cos_omega):
    """
    Equation 24 in Knocke [1].

    Args:
        sin_latitude:    float, sin of the latitude.
        cos_omega:       float, cos of 2*pi*days/365.25, days is number of days since reference epoch Dec 22 1981.
    Returns:
        emissivity:      float, equation 24.
    """
    P0 = 1
    x = sin_latitude
    P1 = x
    P2 = 0.5 * (3 * x**2 - 1)

    return 0.68 * P0 - 0.07 * cos_omega * P1 - 0.18 * P2


def knocke_a(sin_latitude, cos_omega):
    """
    Equation 23 in Knocke [1].

    Args:
        sin_latitude:    float, sin of the latitude.
        cos_omega:       float, cos of 2*pi*days/365.25, days is number of days since reference epoch Dec 22 1981.
    Returns:
        albedo:          float, equation 23.
    """
    P0 = 1
    x = sin_latitude
    P1 = x
    P2 = 0.5 * (3 * x**2 - 1)

    return 0.34 * P0 + 0.10 * cos_omega * P1 - 0.29 * P2
