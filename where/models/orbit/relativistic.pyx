# cython: profile=True
"""Calculates the force acting on the satellite from relativistic effects

Description:

Calculates the force acting on the satellite from relativistic effects, following [1].

References:

   [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS Technical Note No. 36, BKG (2010)

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.math.constant import constant
from midgard.dev import log

# Where imports
from where import apriori

cdef double GM_earth, GM_sun, c, J
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
    global GM_earth
    global GM_sun
    global c, num_param
    global g2i
    global earth_pos
    global earth_vel
    global J 
    
    # Set gravitational constants and speed of light
    GM_earth = constant.get("GM", source="egm_2008")
    GM_sun = constant.get("GM_sun")
    c = constant.get("c")

    # Number of parameters to be estimated
    num_param = len(force_parameters)

    # Position and velocity of the Earth with respect to the Sun
    eph = apriori.get("ephemerides")
    earth_pos = eph.pos_gcrs("earth", time = epochs) - eph.pos_gcrs("sun", time=epochs)
    earth_vel = eph.vel_gcrs("earth", time = epochs) - eph.vel_gcrs("sun", time=epochs)
    
    # Earth's angular momentum per unit mass
    J = constant.get("J")
    
    # Transformation matrices
    g2i = gcrs2itrs


def relativistic(sat_pos_gcrs, sat_vel_gcrs, int current_step, **_not_used):
    """Compute relativistic gravitational force on satellite

    Args:
        sat_pos_gcrs:     Satellite position in GCRS.
        sat_vel_gcrs:     Satellite velocity in GCRS.
        current_step:     Int, step number of current step of integrator.
        _not_used:        Unused variables.

    Returns:
        Acceleration and equation for state transition matrix due to relativistic effects.
    """
    gamma = 1
    beta = 1    
    gcrs2itrs = g2i[current_step]
    r = np.linalg.norm(sat_pos_gcrs)
    v2 = np.dot(sat_vel_gcrs, sat_vel_gcrs) 
    rv = np.dot(sat_pos_gcrs, sat_vel_gcrs)
    rxv = np.cross(sat_pos_gcrs, sat_vel_gcrs)
    # TODO: Not sure if this is the correct treatment of J? 
    J_itrs = np.array([0, 0, J])
    J_gcrs = np.dot(gcrs2itrs.T, J_itrs)
    rJ = np.dot(sat_pos_gcrs, J_gcrs)
    vxJ = np.cross(sat_vel_gcrs, J_gcrs)
    R = earth_pos[current_step, :]
    R_dot = earth_vel[current_step, :]

    # Equation 10.12 from [1]:
    # Schwarzschild terms:
    acc1 = GM_earth / (c**2 * r**3) * ((2 * (beta + gamma) * GM_earth / r - gamma * v2) * sat_pos_gcrs + 2 * (1 + gamma) * rv * sat_vel_gcrs) 
    # Lense-Thirring precession:
    acc2 = (1 + gamma) * GM_earth / (c**2 * r**3) * (3 / r**2 * rxv * rJ + vxJ) 
    # de Sitter precession:
    acc3 = (1 + 2 * gamma) * np.cross(np.cross(R_dot, (-GM_sun * R / (c**2 * np.linalg.norm(R)**3 ))), sat_vel_gcrs)
    # Sum of relativistic effects:
    acc = acc1 + acc2 + acc3

    # Assume negligible state transition matrix
    trans = np.zeros((3, 6))

    # No parameters to be estimated
    sens = np.zeros((3, num_param))
   
    return (acc, trans, sens)
