# cython: profile=True
"""Calculates the empirical force acting on the satellite

Description:

Empirical accelerations, to account for unmodeled forces.

References:

   [1] O. Montenbruck and E. Gill: Satellite Orbits, Springer, 2000, section 3.7.4.


$Revision: 14978 $
$Date: 2018-04-30 19:01:11 +0200 (Mon, 30 Apr 2018) $
$LastChangedBy: hjegei $
"""
# Standard library imports
import sys
import math

# External library imports
import numpy as np

# Where imports
from where import apriori
from where.lib import constant
from where.lib import log
from where.lib import plugins

cdef double GM = constant.GM
cdef double[:] a0, a1, a2

def register_entry_point():
    """Register entry points for setup and later calls."""
    return dict(setup=empirical_setup, call=empirical)


def empirical_setup(rundate, force_parameters, sat_name, time_grid, epochs, body_pos_gcrs, body_pos_itrs, bodies, gcrs2itrs):
    """Set up module variables used later during calculation.
    
    Args: 
        rundate:           Time of integration start
        force_parameters:  Dict of parameters to be estimated
        sat_name:          Name of satellite
        time_grid:         Table of times in seconds since rundate, in utc. 
        epochs:            time_grid converted to Time objects, in utc.
        body_pos:          The positions of the bodies in the solar system, in gcrs
        bodies:            List of bodies
        gcrs2itrs:         List of transformation matrices, one for each time in epochs.
    """
    global GM, a0, a1, a2
    cdef double[:, :] a = np.zeros((3,3))
    cdef int i, j

    for i in range(0, 3):
        for j in range(0, 3):
            key = "a"+str(i)+str(j)
            if key in force_parameters:
                a[i, j] = force_parameters[key]

    a0 = a[0, :]
    a1 = a[1, :]
    a2 = a[2, :]

def empirical(double[:] sat_pos_gcrs, double[:] sat_vel_gcrs, force_parameters, **_not_used):
    """Compute empirical force on satellite

    Args:
        sat_pos_itrs:            Satellite position in GCRS.
        sat_vel_itrs:            Satellite velocity in GCRS.
        sat_name:                Name of satellite.
        force_parameters:        Force parameters to be estimated, dict.
        _not_used:               Unused variables.

    Returns:
        Acceleration and equation for state transition matrix due
        to unmodeled forces.
    """

    cdef double true_anomaly = orbital_elements(sat_pos_gcrs, sat_vel_gcrs)

    cdef double[:] acc
    cdef double[:] pos_unit, vel_unit
    cdef double[:] acc_gcrs
    cdef double[:, :] trans, sens
    cdef double r, v
    cdef int i, j

    # In radial, cross-track and along-track system.
    acc = np.zeros(3)

    for i in range(0, 3):
        acc[i] = a0[i] + a1[i] * math.sin(true_anomaly) + a2[i] * math.cos(true_anomaly)

    cdef double x = sat_pos_gcrs[0]
    cdef double y = sat_pos_gcrs[1]
    cdef double z = sat_pos_gcrs[2]

    r = np.sqrt(x**2 + y**2 + z**2)

    pos_unit = np.zeros(3)
    vel_unit = np.zeros(3)

    for i in range(0, 3):
        pos_unit[i] = sat_pos_gcrs[i] / r

    cdef double xdot = sat_vel_gcrs[0]
    cdef double ydot = sat_vel_gcrs[1]
    cdef double zdot = sat_vel_gcrs[2]

    v = np.sqrt(xdot**2 + ydot**2 + zdot**2)

    for i in range(0, 3):
        vel_unit[i] = sat_vel_gcrs[i] / v

    acc_gcrs = np.zeros(3)

    for i in range(0, 3):
        acc_gcrs[i] = acc[0] * pos_unit[i] + acc[1] * np.cross(vel_unit, pos_unit)[i] + acc[2] * vel_unit[i]

    dacc_dpos = (1 / r**3 * np.matrix([[y**2 + z**2, -x * y, -x * z],
                                       [-x * y, x**2 + z**2, -y * z],
                                       [-x * z, -y * z, x**2 + y**2]]) @
                 (acc[0] * np.eye(3) + acc[1] * np.cross(vel_unit, np.eye(3))))


    dacc_dvel = (1 / v**3 * np.matrix([[ydot**2 + zdot**2, -xdot * ydot, -xdot * zdot],
                                       [-xdot * ydot, xdot**2 + zdot**2, -ydot * zdot],
                                       [-xdot * zdot, -ydot * zdot, xdot**2 + ydot**2]]) @
                 (acc[2] * np.eye(3) + acc[1] * np.cross(np.eye(3), pos_unit)))

    trans = np.hstack((dacc_dpos, dacc_dvel))

    sens = np.zeros((3, len(force_parameters)))
    dacc_dp = np.zeros((3, 9))

    cdef double[:] pos_sin = np.zeros(3)
    cdef double[:] vel_sin = np.zeros(3)
    cdef double[:] cross_sin = np.zeros(3)
    cdef double[:] pos_cos = np.zeros(3)
    cdef double[:] vel_cos = np.zeros(3)
    cdef double[:] cross_cos = np.zeros(3)
    cdef double[:] cross = np.cross(vel_unit, pos_unit)

    for j in range(0, 3):
        pos_sin[j] = pos_unit[j] * math.sin(true_anomaly)
        vel_sin[j] = vel_unit[j] * math.sin(true_anomaly)
        cross_sin[j] = cross[j] * math.sin(true_anomaly)
        pos_cos[j] = pos_unit[j] * math.cos(true_anomaly)
        vel_cos[j] = vel_unit[j] * math.cos(true_anomaly)
        cross_cos[j] = cross[j] * math.cos(true_anomaly)

    # The derivative of the acceleration with respect to the empirical parameters
    dacc_dp[:, 0] = pos_unit
    dacc_dp[:, 1] = cross
    dacc_dp[:, 2] = vel_unit
    dacc_dp[:, 3] = pos_sin
    dacc_dp[:, 4] = cross_sin
    dacc_dp[:, 5] = vel_sin
    dacc_dp[:, 6] = pos_cos
    dacc_dp[:, 7] = cross_cos
    dacc_dp[:, 8] = vel_cos

    # Bookkeeping
    keys = list(force_parameters.keys())
    for i in range(0, len(keys)):
        if keys[i].startswith('a'):
            for j in range(0, 3):
                k1 = int(keys[i][1: 2])
                k2 = int(keys[i][2: 3])
                if k1 == 0:
                    k = k2
                if k1 == 1: 
                    k = k2 + 3
                if k1 == 2:
                    k = k2 + 6
                sens[j, i] = dacc_dp[j, k]

    return (acc_gcrs, trans, sens)


def orbital_elements(double[:] sat_pos, double[:] sat_vel):
    """Computes the true anomaly from satellite position and -velocity

    Args:
        sat_pos:  Satellite position in GCRS.
        sat_vel:  Satellite velocity in GCRS.

    Returns:
        True anomaly:  Angle in radians.

    References: Section 2.2.4 in Montenbruck and Gill [1].
    """
    cdef double[:] h = np.cross(sat_pos, sat_vel)
    cdef double r = np.linalg.norm(sat_pos)
    cdef double sat_vel_squared = 0
    for i in range(0, 3):
        sat_vel_squared += sat_vel[i]**2

    cdef double a = (2 / r - sat_vel_squared / GM)**(-1)

    if a < 0:
        log.fatal("Satellite escaped!!")
        sys.exit(0)

    cdef double p = 0
    for i in range(0, 3):
        p += h[i]**2 / GM

    cdef double e = math.sqrt(1 - p / a)
    cdef double n = math.sqrt(GM / a**3)
    cdef double E = math.atan((np.dot(sat_pos, sat_vel)/(a**2 * n))/(1 - r / a))

    # Make sure we are in the right quadrant:
    if 1 - r / a < 0:
        E += math.pi

    cdef double true_anomaly = math.atan((math.sqrt(1 - e**2) * math.sin(E))/(math.cos(E) - e))

    # Again, make sure we are in the right quadrant:
    if math.cos(E) - e < 0:
        true_anomaly += math.pi

    return true_anomaly
