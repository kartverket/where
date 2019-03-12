# cython: profile=True
"""Calculates the force on the satellite from the gravity field of the Earth

Description:

This model calculates the gravitational force based on gravity field coefficients C and S following
Montenbruck and Gill [1].

References:
[1] Montenbruck, Oliver and Gill, Eberhard, Satellite Orbits, Springer Verlag, 2000.

[2] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS Technical Note No. 36, BKG (2010)

"""

# Standard library imports
import numpy as np
import cython
import math

# Where imports
from where.lib import config
from midgard.math.constant import constant
from where.lib import log
from where import apriori

cdef double GM, R
cdef double[:, :] C, S
cdef int degree_and_order
cdef double[:, :, :] gcrs2itrs
cdef int c20_index


def register_entry_point():
    """Register entry points for setup and later calls."""
    return dict(setup=gravity_earth_setup, call=gravity_earth)


def gravity_earth_setup(
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
        bodies:            The bodies in the solar system
        gcrs2itrs:         List of transformation matrices, one for each time in epochs.
    """
    global C, S
    global GM, R
    global degree_and_order
    global g2i
    global c20_index

    GM = constant.get("GM", source="egm_2008")
    R = constant.get("a", source="egm_2008")
    gravity_field = config.tech.gravity_field.str
    truncation_level = config.tech.gravity_truncation_level.int
    gravity_coeffs = apriori.get("gravity", gravity_field=gravity_field, truncation_level=truncation_level,
                                 rundate=rundate)
    C = gravity_coeffs["C"]
    S = gravity_coeffs["S"]
    #Assume for now that degree and order of gravity field are equal.
    degree_and_order = C.shape[0]
    g2i = gcrs2itrs
    c20_index = -1
    if "c20" in force_parameters:
        C[2, 0] = force_parameters["c20"]
        c20_index = list(force_parameters.keys()).index("c20")


def gravity_earth(double[:] sat_pos_itrs, force_parameters, int current_step, **_not_used):
    """Compute force on satellite from the gravity field of the earth

    Following section 3.2.5 in Montenbruck and Gill [1] we calculate the acceleration caused by the gravity field of
    the earth from the gravity field coefficients \f$ C_{nm} \f$ and \f$ S_{nm} \f$, in addition to the Legendre terms
    \f$ V_{nm} \f$ and \f$ W_{nm} \f$ using equation (3.33).

    In addition, we calculate the transition matrix \f$ \phi \f$ using equation (7.42),

    \f[ \dot{\vec \phi}(t, t_0) = \left( \begin{array}{cc} 0_{3 \times 3} &
    I_{3 \times 3} \\ \frac{\partial \ddot{\vec r}(\vec r, \dot{\vec r}, t)}{
    \partial \vec r(t)} & \frac{\partial \ddot{\vec r}(\vec r, \dot{\vec r},
    t)}{\partial \dot{\vec r}(t)} \end{array} \right) \cdot \vec \phi(t, t_0)
    . \f]

    The acceleration with respect to position (lower left corner) is calculated using the equations (7.65) to (7.69).

    The terms \f$ V, W \f$ are required up to degree and order 2 higher than the truncation level of the gravity field,
    as can be seen from the equations for the transition matrix.

    The calculation of the gravity field of the earth is done in the earth fixed ITRS system. the coordinates are
    transformed to the space-fixed GCRS system before they are returned, according to equations (3.34) and (7.70).

    Args:
        sat_pos_itrs:      Satellite position in ITRS.
        force_parameters:  Force parameters to be estimated.
        current_step:      Int, step number of current step of integrator.
        _not_used:         Unused variables.

    Returns:
        Acceleration and transition matrix due to earth gravity field in GCRS.
    """
    cdef int i, f, n, m
    cdef double dxx = 0, dxy = 0, dxz = 0, dyz = 0, dzz = 0
    cdef double[:] acc = np.zeros(3)
    cdef double[:, :] V, W, trans, trans_gcrs
    cdef double[:, :] sens_itrs, sens_gcrs
    cdef double[:, :, :] VW
    VW = compute_VW(sat_pos_itrs)

    V = VW[:, :, 0]
    W = VW[:, :, 1]

    # Acceleration forces, equation (3.33)
    # Transition matrix, equations (7.65) - (7.70)

    for n in range(0, degree_and_order):
        for m in range(0, n + 1):
            f = (n - m + 2) * (n - m + 1)    # Scaling factor
            # Denormalize the C and S coefficients in the equations
            normalization_factor = np.sqrt(math.factorial(n - m) * (2 * n + 1) * (2 - (m == 0)) /
                                           math.factorial(n + m))
            # The m = 0 case is handled separately:
            if m == 0:
                acc[0] += -normalization_factor * C[n, 0] * V[n+1, 1]
                acc[1] += -normalization_factor * C[n, 0] * W[n+1, 1]
                acc[2] += normalization_factor * (n + 1) * (-C[n, 0] * V[n+1, 0] - S[n, 0] * W[n+1, 0])
                dxx += normalization_factor * (C[n, 0] * V[n+2, 2] - f * C[n, 0] * V[n+2, 0]) / 2
                dxy += normalization_factor * C[n, 0] * W[n+2, 2] / 2
                dxz += normalization_factor * (n + 1) * C[n, 0] * V[n+2, 1]
                dyz += normalization_factor * (n + 1) * C[n, 0] * W[n+2, 1]
                dzz += normalization_factor * f * (C[n, 0] * V[n+2, 0] + S[n, 0] * W[n+2, 0])
                continue

            # For some derivatives also m=1 needs special treatment
            if m == 1:
                dxx += normalization_factor * (C[n, 1] * V[n + 2, 3] + S[n, 1] * W[n + 2, 3]
                        + f * (-3 * C[n, 1] * V[n + 2, 1] - S[n, 1] * W[n+2, 1])) / 4
                dxy += normalization_factor * (C[n, 1] * W[n + 2, 3] - S[n, 1] * V[n + 2, 3]
                        + f * (-C[n, 1] * W[n + 2, 1] - S[n, 1] * V[n + 2, 1])) / 4
            else:
                dxx += (normalization_factor
                        * (C[n, m] * V[n + 2, m + 2] + S[n, m] * W[n + 2, m + 2] + 2 * (n - m + 2) * (n - m + 1)
                           * (-C[n, m] * V[n + 2, m] - S[n, m] * W[n + 2, m]) + (n - m + 4) * (n - m + 3) * (n - m + 2)
                           * (n - m + 1) * (C[n, m] * V[n + 2, m - 2] + S[n, m] * W[n + 2, m - 2])) / 4)
                dxy += (normalization_factor
                        * (C[n, m] * W[n + 2, m + 2] - S[n, m] * V[n + 2, m + 2] + (n - m + 4) * (n - m + 3)
                           * (n - m + 2) * (n - m + 1) * (-C[n, m] * W[n + 2, m - 2] + S[n, m] * V[n + 2, m - 2])) / 4)

            # For 0 < m <= n:
            f = (n - m + 2) * (n - m + 1)    # Scaling factor
            acc[0] += normalization_factor * ((-C[n, m] * V[n + 1, m + 1] - S[n, m] * W[n + 1, m + 1])
                       + f * (C[n, m] * V[n + 1, m - 1] + S[n, m] * W[n + 1, m - 1])) / 2
            acc[1] += normalization_factor * ((-C[n, m] * W[n + 1, m + 1] + S[n, m] * V[n + 1, m + 1])
                       + f * (-C[n, m] * W[n + 1, m - 1] + S[n, m] * V[n + 1, m - 1])) / 2
            acc[2] += normalization_factor * (n - m + 1) * (-C[n, m] * V[n + 1, m] - S[n, m] * W[n + 1, m])

            dxz += (normalization_factor
                    * ((n - m + 1) * (C[n, m] * V[n + 2, m + 1] + S[n, m] * W[n + 2, m + 1]) + (n - m + 3) * f
                       * (-C[n, m] * V[n + 2, m - 1] - S[n, m] * W[n + 2, m - 1])) / 2)
            dyz += (normalization_factor
                    * ((n - m + 1) * (C[n, m] * W[n + 2, m + 1] - S[n, m] * V[n + 2, m + 1]) + (n - m + 3) * f
                       * (C[n, m] * W[n + 2, m - 1] - S[n, m] * V[n + 2, m - 1])) / 2)
            dzz += normalization_factor * f * (C[n, m] * V[n + 2, m] + S[n, m] * W[n + 2, m])

    for i in range(3):
        acc[i] *= GM / R**2

    cdef double fact = GM / R**3

    dxx *= fact
    dxy *= fact
    dxz *= fact
    dzz *= fact
    dyz *= fact

    trans = np.array([[dxx, dxy, dxz], [dxy, -dxx - dzz, dyz], [dxz, dyz, dzz]])
    trans_gcrs = np.zeros((3, 6))
    gcrs2itrs = g2i[current_step]

    # In numpy:
    # trans_gcrs = np.hstack((np.dot(gcrs2itrs.T, np.dot(trans, gcrs2itrs)), np.zeros((3, 3))))
    for i in range(0, 3):
        for j in range(0, 3):
            for k in range(0, 3):
                for l in range(0, 3):
                    trans_gcrs[i, j] += gcrs2itrs[k, i] * trans[k, l] * gcrs2itrs[l, j]

    sens_itrs = np.zeros((3, len(force_parameters)))
    sens_gcrs = np.zeros((3, len(force_parameters)))
    if not c20_index == -1:
        # Equation (7.73) in [1].
        sens_itrs[0, c20_index] = -GM / R**2 * V[3, 1]
        sens_itrs[1, c20_index] = -GM / R**2 * W[3, 1]
        sens_itrs[2, c20_index] = -GM / R**2 * V[3, 0]
    for i in range(0, 3):
        for j in range(0, len(force_parameters)):
            for k in range(0, 3):
                sens_gcrs[i, j] += gcrs2itrs[k, i] * sens_itrs[k, l]

    # Transform to space fixed system before returning, eqs (3.34) and (7.70)
    return (np.dot(gcrs2itrs.T, acc), trans_gcrs, sens_gcrs)


cdef double[:, :, :] compute_VW(double[:] pos_xyz):
    """Computing the V- and W-coefficients V and W recursively

    The V- and W-coefficients are based on Legendre polynomials, and used when calculating the gravity potential. The
    coefficients are calculated using recurrence relations as described in section 3.2.4 of Montenbruck and Gill [1].

    Args:
        pos_xyz:          Position of satellite as a 3-vector.

    Returns:
        Two matrices V and W with coefficients.
    """
    cdef double[:, :, :] VW
    cdef double x, y, z, r, f
    cdef int n, m
    VW = np.zeros((degree_and_order + 2, degree_and_order + 2, 2))

    x = pos_xyz[0]
    y = pos_xyz[1]
    z = pos_xyz[2]
    r = (x**2 + y**2 + z**2)**0.5

    if r < R:
        log.fatal("SATELLITE CRASHED !!!")
    elif r > 36e6:
        log.warn(f"Satellite flying high, r={r} m")
        if r > 40e7:
            log.fatal("SATELLITE FLYING TOO HIGH, BYE BYE SATELLITE")

    f = R / r**2   # Common factor

    VW[0, 0, 0] = R / r
    VW[1, 0, 0] = z * f * VW[0, 0, 0]

    # First compute the zonal terms V[n,0]. The terms W[n,0] are always zero.
    for n in range(2, degree_and_order + 2):
        VW[n, 0, 0] = f * ((2 * n - 1) * z * VW[n - 1, 0, 0] - (n - 1) * R * VW[n - 2, 0, 0]) / n

    for m in range(1, degree_and_order + 2):
        # Compute the diagonal matrix elements, called the tesseral terms.
        VW[m, m, 0] = (2 * m - 1) * f * (x * VW[m - 1, m - 1, 0] - y * VW[m - 1, m - 1, 1])
        VW[m, m, 1] = (2 * m - 1) * f * (x * VW[m - 1, m - 1, 1] + y * VW[m - 1, m - 1, 0])

        # Compute the remaining terms
        for n in range(m + 1, degree_and_order + 2):
            VW[n, m, 0] = f * ((2 * n - 1) * z * VW[n - 1, m, 0] - (n + m - 1) * R * VW[n - 2, m, 0]) / (n - m)
            VW[n, m, 1] = f * ((2 * n - 1) * z * VW[n - 1, m, 1] - (n + m - 1) * R * VW[n - 2, m, 1]) / (n - m)

    return VW
