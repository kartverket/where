"""Calculate the VLBI geometric delay for near field targets

Description:
------------

Calculate the geometric delay using a near field model as described in :cite:`jaron2019` and :cite:`deuv2012`.

The model derived in `jaron2019` is intended for Earth satellites and uses a linear approximation for short 
term station and satellite motion. This approximation allows for a analytical solution and the equations are 
expressed in the GCRS. The gravitational effect of celestial bodies on the delay is described in `deuv2012`.

"""

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.constant import constant

# Where imports
from where import apriori
from where.data.position import PosVel
from where.data.time import TimeDelta



@plugins.register
def vlbi_near_field_delay(dset):
    r"""Calculate the theoretical delay dependent on the baseline

    TODOTODOTODOTODO :
    --------------------------------------
    The implementation is described in IERS Conventions :cite:`iers2010`, section 11.1, in particular equation
    (11.9). We do not take the gravitational delay into account here (see
    :mod:`where.models.delay.vlbi_gravitational_delay`), and multiply by :math:`c` to get the correction in
    meters. Thus, we implement the following equation:

    .. math::
       \mathrm{correction} = \frac{- \hat K \cdot \vec b \bigl[ 1 - \frac{(1 + \gamma) U}{c^2}
                             - \frac{| \vec V_\oplus |^2}{2 c^2} - \frac{\vec V_\oplus \cdot \vec w_2}{c^2} \bigr]
                             - \frac{\vec V_\oplus \cdot \vec b}{c} \bigl[ 1
                             + \frac{\hat K \cdot \vec V_\oplus}{2 c} \bigr]}{1
                             + \frac{\hat K \cdot (\vec V_\oplus + \vec w_2)}{c}}

    with

    * :math:`\hat K` -- the unit vector from the barycenter to the source in the absence of gravitational or
      aberrational bending,

    * :math:`\vec b` -- the GCRS baseline vector at the time :math:`t_1` of arrival, :math:`\vec x_2(t_1) - \vec
      x_1(t_1)`,

    * :math:`\gamma` -- the parameterized post-Newtonian (PPN) gamma, equal to 1 in general relativity theory,

    * :math:`U` -- the gravitational potential at the geocenter, neglecting the effects of the Earth's mass. At the
      picosecond level, only the solar potential need be in included in :math:`U` so that :math:`U = G M_\odot / | \vec
      R_{\oplus_\odot} |` where :math:`\vec R_{\oplus_\odot}` is the vector from the Sun to the geocenter,

    * :math:`\vec V_\oplus` -- the barycentric velocity of the geocenter,

    * :math:`\vec w_2` -- the geocentric velocity of station 2.

    Each term in the correction is calculated in separate functions. and stored in the Dataset in a table called
    ``vlbi_vacuum_delay``.
    -------------------------------------------------------

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array: Projected baseline in meters for each observation.

    """
    orbit = apriori.get("simple_orbit", rundate=dset.analysis["rundate"], days_before=0, days_after=1)
    # TODO
    # idx_sat = True when observation is to a satellite

    delay = np.zeros(dset.num_obs)
    
    # Apriori values given at epoch t1
    t1 = dset.time.tcg
    x1_t1 = dset.site_pos_1.gcrs.pos # station_1 at epoch t1
    x2_t1 = dset.site_pos_2.gcrs.pos # station_2 at epoch t1
    x0_t1 = dset.sat_pos.gcrs.pos # satellite position at epoch t1
    v0_t1 = dset.sat_pos.gcrs.vel.val # satellite velocity at at epoch t1
    v2_t1 = dset.site_pos_2.gcrs.vel.val # station_2 velocity at epoch t1
    
    # First approximation to light travel time
    delta1 = (x1_t1 - x0_t1).length / constant.c # eq. 4 # seconds
    delta2 = (x2_t1 - x0_t1).length / constant.c # eq. 4 # seconds
    
    # Convert to TimeDelta objects
    delta1 = TimeDelta(delta1, fmt="seconds", scale="tcg")
    delta2 = TimeDelta(delta2, fmt="seconds", scale="tcg")
    
    t0_tilde = t1 - delta1 # approximation to t0
    tau_tilde = delta2 - delta1 # eq. 7
    t2_tilde = t1 + tau_tilde # approximation to t2
 
    # TODO "loop over satellites in session
    sat_pos = orbit["G10"]["pos"](t0_tilde)
    sat_vel = orbit["G10"]["vel"](t0_tilde)
    sat_posvel = PosVel(np.concatenate((sat_pos, sat_vel), axis=1), system="trs", time=t0_tilde)
    
    # Satellite position and velocity at t0
    x0_t0_tilde = sat_posvel.gcrs.pos.val
    v0_t0_tilde = sat_posvel.gcrs.vel.val
    
    # Linearized satellite position at t1
    dt_10 = (t1 - t0_tilde).seconds[:, None]
    x0_bar_t1 = v0_t0_tilde * dt_10 + x0_t0_tilde # eq. 5
    
    gamma0_2 = np.sqrt(1 - (v0_t1[:, None, :] @ v0_t1[:, :, None])[:, 0, 0] / constant.c**2) # eq. 15
    x01 = x0_bar_t1 - x1_t1.val # eq. 16
    
    # TODO
    t_g01 = 0
    
    # eq. 14
    x_dot_v_1 = (x01[:, None, :] @ v0_t1[:, :, None])[:, 0, 0] / constant.c ** 2 # Intermediate variable
    x01_dot_x01 = (x01[:, None, :] @ x01[:, :, None])[:, 0, 0] # Intermediate variable
    # Time of emmison of the signal relative to t1
    delta_t0 = gamma0_2 * (x_dot_v_1 - t_g01) - \
        np.sqrt(gamma0_2 ** 2 * (x_dot_v_1 - t_g01) ** 2 + gamma0_2 * (x01_dot_x01 / constant.c ** 2 - t_g01 ** 2))
    
    
    # Assume station_2 has no motion beweteen t1 and t2 in a terrestrial reference system
    site_pos_2_t2 = PosVel(dset.site_pos_2.val, system="trs", time=t2_tilde)
    # GCRS posiion at t2
    x2_t2_tilde = site_pos_2_t2.gcrs.pos.val 
    v2_t2_tilde = site_pos_2_t2.gcrs.vel.val
    
    # Linearized station_2 position at t1
    dt_12 = (t1 - t2_tilde).seconds[:, None]
    x2_bar_t1 = v2_t2_tilde * dt_12 + x2_t2_tilde # eq.8
    
    gamma2_2 = np.sqrt(1 - (v2_t1[:, None, :] @ v2_t1[:, :, None])[:, 0, 0] / constant.c**2) # eq. 18 
    x02 = x0_bar_t1 - x2_bar_t1 + (v0_t1 - v2_t1) * delta_t0[:, None] # eq. 19
    
    # TODO
    t_g02 = 0
    
    # eq. 17
    x_dot_v_2 = (x02[:, None, :] @ v2_t1[:, :, None])[:, 0, 0] / constant.c ** 2 # Intermediate variable
    x02_dot_x02 = (x02[:, None, :] @ x02[:, :, None])[:, 0, 0] # Intermediate variable
    # Time of reception of the signal relative to t1
    delta_t2 = - gamma2_2 * (x_dot_v_2 - t_g02) + \
        np.sqrt(gamma2_2 ** 2 * (x_dot_v_2 - t_g02) ** 2 + gamma2_2 * (x02_dot_x02 / constant.c ** 2 - t_g02 ** 2))
    
    # Convert from TCG to TT
    delay = (delta_t2 + delta_t0) * (1 - constant.L_G) # eq. 10 

    return delay * constant.c


def grav_delay(dset):
    """
    Eq. 14 in :cite:`deuv2012`
    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array: Gravitational delay in meters for each observation.
    """
    
    eph = apriori.get("ephemerides", time=dset.time)
    grav_delay = np.zeros(dset.num_obs)
    
    # List of celestial bodies. Major moons are also recommended, like Titan, Ganymedes, ...
    bodies = [
        "mercury barycenter",
        "venus barycenter",
        "earth",
        "moon",
        "mars barycenter",
        "jupiter barycenter",
        "saturn barycenter",
        "uranus barycenter",
        "neptune barycenter",
        "pluto barycenter",
        "sun",
    ]
    
    #TODO
    return 0
    
    gamma = 1 # PPN parameter. Equal to 1 in general relativity
    
    factor = (1 + gamma)/constant.c**2 # Multiplication factor repeated frequently in equation
    GM_sun = constant.get("sun", source=eph.ephemerides)
    
    S_factor = factor * constant.G * GM_sun # S is short for Sun
    # index 0 is satelitte
    # index 1 or 2 is station
    # index 01 and 02 is vector between station and satellite
    # Capital R is barycentric position/vector
    
    # TODO: Confirm this. What about time scale?
    # Ref. IERS 2010 conventions eq. 11.6 used for grav delay in consensus model
    # Transformation between BCRS and GRCS is approximated by X_bcrs = X_earth + X_gcrs
    R_S_T0 = eph.pos_bcrs("sun")
    R_S_T1 = eph.pos_brcs("sun")
    R_1 = eph.pos_bcrs("earth") + dset.site_pos_1.gcrs.pos.val
    R_0 = eph.pos_bcrs("earth") + dset.sat_pos.gcrs.pos.val
    
    R_01 = R_1 - R_0
    R_0_S = R_0 - R_S # eq. 16
    R_1_S = R_1 - R_S # eq. 16
    grav_delay_sun = S_factor/constant.c * \
        np.log((R_0_S + R_1_S + R_01_S + S_factor)/(R_0_S + R_1_S - R_01_S + S_factor))

    grav_delay_bodies = 0

# def term_1(dset, proj_Kb, _, _ve):
#     r"""Main part of the vacuum delay is the baseline in the source direction
#
#     The term :math:`\hat K \cdot \vec b \cdot \bigl( -1 \bigr)` scaled by the denominator :math:`1 + \frac{\hat K \cdot
#     (\vec V_\oplus + \vec w_2)}{c}`.
#
#     Args:
#         dset:    Model input data.
#         proj_Kb: Scaled projection of baseline in direction of source.
#
#     Returns:
#         Numpy array: Part of vacuum delay.
#     """
#     return -proj_Kb
#
#
# def term_2(dset, proj_Kb, _, _ve):
#     r"""Part of the vacuum delay dependent on the gravitational potential
#
#     The term :math:`\hat K \cdot \vec b \cdot \frac{(1 + \gamma) U}{c^2}` scaled by the denominator :math:`1 +
#     \frac{\hat K \cdot (\vec V_\oplus + \vec w_2)}{c}`.
#
#     The parameterized post-Newtonian (PPN) gamma, :math:`\gamma` is equal to 1 in general relativity theory.
#
#     The gravitational potential at the geocenter, \f$ U \f$, neglecting the effects of the Earth's mass is
#     calculated. Following table 11.1 in IERS Conventions [2], only the solar potential need to be included at the
#     picosecond level. That is
#
#     \f[ U = G M_\odot / | \vec R_{\oplus_\odot} | \f]
#
#     where \f$ \vec R_{\oplus_\odot} \f$ is the vector from the Sun to the geocenter. We calculate the latter using the
#     ephemerides.
#
#     Args:
#         dset:    Model input data.
#         proj_Kb: Scaled projection of baseline in direction of source.
#
#     Returns:
#         Numpy array: Part of vacuum delay.
#     """
#     gamma = 1.0
#     eph = apriori.get("ephemerides", time=dset.time)
#     grav_potential = constant.GM_sun / np.linalg.norm(eph.pos_gcrs("sun"), axis=1)
#     return proj_Kb * (1 + gamma) * grav_potential / constant.c ** 2
#
#
# def term_3(dset, proj_Kb, _, vel_earth):
#     r"""Correction to delay based on earth's movement in space
#
#     The term :math:`\hat K \cdot \vec b \cdot \frac{| \vec V_\oplus |^2}{2 c^2}` scaled by the denominator :math:`1 +
#     \frac{\hat K \cdot (\vec V_\oplus + \vec w_2)}{c}`.
#
#     Args:
#         dset:    Model input data.
#         proj_Kb: Scaled projection of baseline in direction of source.
#
#     Returns:
#         Numpy array: Part of vacuum delay.
#     """
#     return proj_Kb * 0.5 * (vel_earth[:, None, :] @ vel_earth[:, :, None] / constant.c ** 2)[:, 0, 0]
#
#
# def term_4(dset, proj_Kb, _, vel_earth):
#     r"""Correction to the delay caused by earth's rotation
#
#     The term :math:`\hat K \cdot \vec b \cdot \frac{\vec V_\oplus \cdot \vec w_2}{c^2}` scaled by the denominator
#     :math:`1 + \frac{\hat K \cdot (\vec V_\oplus + \vec w_2)}{c}`.
#
#     Args:
#         dset:    Model input data.
#         proj_Kb: Scaled projection of baseline in direction of source.
#
#     Returns:
#         Numpy array: Part of vacuum delay.
#     """
#     return proj_Kb * (vel_earth[:, None, :] @ dset.site_pos_2.gcrs.vel.mat / constant.c ** 2)[:, 0, 0]
#
#
# def term_5(dset, _, proj_Vb, _ve):
#     r"""Part of the delay due to earth's movement in space
#
#     The term :math:`- \frac{\vec V_\oplus \cdot \vec b}{c} \cdot \bigl( 1 \bigr)` scaled by the denominator :math:`1 +
#     \frac{\hat K \cdot (\vec V_\oplus + \vec w_2)}{c}`.
#
#     Args:
#         dset:    Model input data.
#         proj_Vb: Scaled projection of baseline in direction of earth's movement.
#
#     Returns:
#         Numpy array: Part of vacuum delay.
#     """
#     return -proj_Vb
#
#
# def term_6(dset, _, proj_Vb, vel_earth):
#     r"""Correction to earth's movement in space
#
#     The term :math:`- \frac{\vec V_\oplus \cdot \vec b}{c} \cdot \frac{\hat K \cdot \vec V_\oplus}{2 c}` scaled by the
#     denominator :math:`1 + \frac{\hat K \cdot (\vec V_\oplus + \vec w_2)}{c}`.
#
#     Args:
#         dset:    Model input data.
#         proj_Vb: Scaled projection of baseline in direction of earth's movement.
#
#     Returns:
#         Numpy array: Part of vacuum delay.
#     """
#     return -proj_Vb * 0.5 * (dset.src_dir.unit_vector[:, None, :] @ vel_earth[:, :, None] / constant.c)[:, 0, 0]
