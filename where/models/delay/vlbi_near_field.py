"""Calculate the VLBI geometric delay for near field targets

Description:
------------

Calculate the geometric delay using a near field model as described in :cite:`jaron2019` and :cite:`deuv2012`.

The model derived in `jaron2019` is intended for Earth satellites and uses a linear approximation for short 
term station and satellite motion. This approximation allows for a analytical solution and the equations are 
expressed in the GCRS. The gravitational effect of celestial bodies on the delay is described in `deuv2012`.

The model derived in `deuv2012` is expressed in the barycentric reference frame and is based on an iterative
solution of the light time equations. This model is valid for the entire solar system. 

"""
# Standard library imports
from datetime import datetime
from typing import Any, Callable, Dict, List, Tuple

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.constant import constant
from midgard.math.unit import Unit

# Where imports
from where import apriori
from where.data.position import PosVel
from where.data.time import Time, TimeDelta
from where.lib import config
from where.lib import log

# Constants for shorter equations
GAMMA = 1 # PPN parameter. Equal to 1 in general relativity 
C = constant.c
L_G = constant.L_G
L_C = constant.L_C

# The name of the satellite in the orbit file is not the same as the name in the NGS testfiles
# Create a small translation table
ngs_to_sp3 = dict()
ngs_to_sp3["GEN-01"] = "L01"
ngs_to_sp3["LAGEOS-1"] = "L51"
ngs_to_sp3["SENTI-6A"] = "L40"

MODEL = __name__.split(".")[-1]

MODELS = {}

def register_model(model: Callable) -> Callable:
    MODELS[model.__name__] = model
    return model

@plugins.register
def vlbi_near_field(dset):
    r"""Calculate the theoretical delay dependent on the baseline

    -------------------------------------------------------

    Args:
        dset:     A Dataset containing model data.

    Returns:

        Numpy array: Near field delay for each observation

    """
    near_field_model = config.tech[MODEL].model.str
    #j_delay = jaron2019(dset)
    #d_delay = deuv2012(dset)
    #import IPython; IPython.embed()
    if near_field_model in MODELS:
        return MODELS[near_field_model](dset)
    else:
        log.error(f"Unknown model {near_field_model} for {MODEL} delay model")

@register_model    
def jaron2019(dset):
    """ 
    Implementation based on
        Jaron et al. 2019: Modelling the VLBI delay for Earth satellites

    This model is expressed in the GCRS and uses TGC-compatible values. This has some consequences (see Hakan et al. 2026
    for explainations):
        Gravitational effects on the delay should only account for the planet Earth
        Apriori values for station and satellite positions are usually TT-compatible and needs to be scaled properly
        The final output needs to be scaled to TT
    """
    # This model is only applicable for near field observations
    idx = dset.near_field_obs
    num_sat_obs = np.sum(idx)
    if num_sat_obs == 0:
        # Skip this model if there are no near field observations
        return np.zeros(dset.num_obs)

    # Apriori values given at epoch t1 (time of arrival for signal at station 1)
    time = dset.time[idx]
    t1 = time.tcg
    x1_t1 = dset.site_pos_1.gcrs.pos[idx] # station_1 position at epoch t1
    x2_t1 = dset.site_pos_2.gcrs.pos[idx] # station_2 position at epoch t1
    x0_t1 = dset.sat_pos.gcrs.pos[idx] # satellite position at epoch t1
    v0_t1 = dset.sat_pos.gcrs.vel.val[idx] / ((1 - L_G)) # satellite velocity at at epoch t1
    v2_t1 = dset.site_pos_2.gcrs.vel.val[idx] / ((1 - L_G)) # station_2 velocity at epoch t1

    # First approximation to light travel time
    delta1 = (x1_t1 - x0_t1).length / C / ((1 - L_G)) # eq. 4 # seconds
    delta2 = (x2_t1 - x0_t1).length / C / ((1 - L_G)) # eq. 6 # seconds
    
    # Convert to TimeDelta objects
    delta1 = TimeDelta(delta1, fmt="seconds", scale="tcg")
    delta2 = TimeDelta(delta2, fmt="seconds", scale="tcg")

    t0_tilde = t1 - delta1 # approximation to t0 (time of emission of signal from satellite)
    tau_tilde = delta2 - delta1 # eq. 7
    t2_tilde = t1 + tau_tilde # approximation to t2 (time of arrival for signal at station 2)

    sat_posvel = _sat_posvel(dset, t0_tilde) 
    
    # Satellite position and velocity at t0
    x0_t0_tilde = sat_posvel.gcrs.pos.val / ((1 - L_G))
    v0_t0_tilde = sat_posvel.gcrs.vel.val / ((1 - L_G))
    
    # Linearized satellite position at t1
    dt_10 = delta1.seconds[:, None]
    x0_bar_t1 = v0_t0_tilde * dt_10 + x0_t0_tilde # eq. 5
    
    gamma0_2 = 1/(1 - (v0_t1[:, None, :] @ v0_t1[:, :, None])[:, 0, 0] / C ** 2) # eq. 15
    x01 = x0_bar_t1 - x1_t1.val # eq. 16

    # Compute t_g01: Relativistic effects on delay from satellite to station 1
    # Based on Deuv, et al (2012) eq. 14, 16, 17
    # Equations are in BCRS. Ephemerides use TDB.
    bodies = ["earth"] # This model is expressed in GCRS and should only include the gravitational effect from the Earth
    R0_T0 = _g2b_pos(x0_t0_tilde, t0_tilde) # satellite position at t0 in BCRS
    R1_T1 = _g2b_pos(x1_t1.val, t1) # station_1 pos at t1 in BCRS
    t_g01_TDB = _deuv_relativistic_term(R0_T0, R1_T1, t0_tilde, time, bodies)
    t_g01 = t_g01_TDB / (1 - L_G)

    # Save TT(=TDB) value to dset
    _save_float_to_dset(dset, idx, f"{MODEL}.grav_1", t_g01_TDB * C, unit="meter", write_level="detail")

    # eq. 14 in jaron2019
    x01_dot_v0 = (x01[:, None, :] @ v0_t1[:, :, None])[:, 0, 0] / C ** 2 # Intermediate variable
    x01_dot_x01 = (x01[:, None, :] @ x01[:, :, None])[:, 0, 0]  / C ** 2 # Intermediate variable
    # Time of emmison of the signal relative to t1
    delta_t0 = gamma0_2 * (x01_dot_v0 - t_g01) - \
        np.sqrt(gamma0_2 ** 2 * (x01_dot_v0 - t_g01) ** 2 + gamma0_2 * (x01_dot_x01 - t_g01 ** 2))
    
    
    # Assume station_2 has no motion beweteen t1 and t2 in a terrestrial reference system
    site_pos_2_t2 = PosVel(dset.site_pos_2.val[idx], system="trs", time=t2_tilde)
    # GCRS posiion at t2
    x2_t2_tilde = site_pos_2_t2.gcrs.pos.val / ((1 - L_G))
    v2_t2_tilde = site_pos_2_t2.gcrs.vel.val / ((1 - L_G))
    
    # Linearized station_2 position at t1
    dt_12 = - tau_tilde.seconds[:, None]
    x2_bar_t1 = v2_t2_tilde * dt_12 + x2_t2_tilde # eq.8
    
    gamma2_2 = 1/(1 - (v2_t1[:, None, :] @ v2_t1[:, :, None])[:, 0, 0] / C ** 2) # eq. 18 
    x02 = x0_bar_t1 - x2_bar_t1 + (v0_t1 - v2_t1) * delta_t0[:, None] # eq. 19
    
    # Compute t_g02: Relativistic effects on delay from satellite to station 2
    # Based on Deuv, et al (2012) eq. 14, 16, 17
    # Equations are in BCRS. Ephemerides use TDB.  
    R2_T2 = _g2b_pos(x2_t2_tilde, t2_tilde) # station_2 pos at t2 in BCRS
    t_g02_TDB = _deuv_relativistic_term(R0_T0, R2_T2, t2_tilde, time, bodies)
    t_g02 = t_g02_TDB / (1 - L_G)
    
    # Save TT(=TDB) value to dset  
    _save_float_to_dset(dset, idx, f"{MODEL}.grav_2", t_g02_TDB * C, unit="meter", write_level="detail")
    
    # eq. 17 in jaron2019
    x02_dot_v2 = (x02[:, None, :] @ v2_t1[:, :, None])[:, 0, 0] / C ** 2 # Intermediate variable
    x02_dot_x02 = (x02[:, None, :] @ x02[:, :, None])[:, 0, 0] / C ** 2 # Intermediate variable

    # Time of reception of the signal relative to t1
    delta_t2 = - gamma2_2 * (x02_dot_v2 - t_g02) + \
        np.sqrt(gamma2_2 ** 2 * (x02_dot_v2 - t_g02) ** 2 + gamma2_2 * (x02_dot_x02 - t_g02 ** 2))
    
    # Convert from TCG to TT
    delay = (delta_t2 + delta_t0) * (1 - L_G) # eq. 10 


    # Save intermediate variables to dataset for reuse in computation of partials
    # All variables are TCG compatible
    _save_float_to_dset(dset, idx, f"{MODEL}.v0_t0_tilde", v0_t0_tilde, unit="(m/s, m/s, m/s)")
    _save_float_to_dset(dset, idx, f"{MODEL}.v2_t2_tilde", v2_t2_tilde, unit="(m/s, m/s, m/s)")
    _save_float_to_dset(dset, idx, f"{MODEL}.x01", x01, unit="(m, m, m)")
    _save_float_to_dset(dset, idx, f"{MODEL}.x02", x02, unit="(m, m, m)")
    _save_float_to_dset(dset, idx, f"{MODEL}.gamma0", np.sqrt(gamma0_2), unit="dimensionless")
    _save_float_to_dset(dset, idx, f"{MODEL}.gamma2", np.sqrt(gamma2_2), unit="dimensionless")
    _save_float_to_dset(dset, idx, f"{MODEL}.delta_t0", delta_t0, unit="seconds")
    # Add more variables for debugging purposes
    _save_time_to_dset(dset, idx, f"{MODEL}.t0_tilde", t0_tilde, write_level="detail")
    _save_time_to_dset(dset, idx, f"{MODEL}.t2_tilde", t2_tilde, write_level="detail")
    _save_float_to_dset(dset, idx, f"{MODEL}.delta_t2", delta_t2, write_level="detail", unit="seconds")
    _save_posvel_to_dset(dset, idx, f"{MODEL}.j_sat_posvel", sat_posvel, write_level="detail")
    
    ## For debugging. See if satellite is above horizon for both stations
    debug = False
    if debug:
        e1 = dset.site_pos_1.elevation_to(dset.sat_pos)
        e2 = dset.site_pos_2.elevation_to(dset.sat_pos)
        sat_visible = (e1 > 0) & (e2 > 0)
        
        dset.add_bool("sat_visible", sat_visible)
    
        import matplotlib.pyplot as plt;
        for bl in dset.unique("baseline"):
            bl_idx = dset.filter(baseline=bl)
            alpha = np.ones(np.sum(bl_idx))
            alpha[dset.sat_visible[bl_idx] == False] = 0.1
            for body in bodies + ["sun"]:
                plt.scatter(dset.time.datetime[bl_idx], dset[f"{MODEL}.grav_{body}_1"][bl_idx]/C, alpha=alpha, label=f"{body}_1")
                plt.scatter(dset.time.datetime[bl_idx], dset[f"{MODEL}.grav_{body}_2"][bl_idx]/C, alpha=alpha, label=f"{body}_2")
            plt.legend(ncol=2, loc='center left', bbox_to_anchor=(1, 0.5))
            plt.title(bl)
            plt.tight_layout()
            plt.show()
            
            for body in bodies + ["sun"]:
                y = (dset[f"{MODEL}.grav_{body}_2"][bl_idx] - dset[f"{MODEL}.grav_{body}_1"][bl_idx])/C
                plt.scatter(dset.time.datetime[bl_idx], y, alpha=alpha, label=f"diff_{body}")
            plt.legend(ncol=1, loc='center left', bbox_to_anchor=(1, 0.5))
            plt.title(bl)
            plt.tight_layout()
            plt.show()

    
        #import IPython; IPython.embed()
    output = np.zeros(dset.num_obs)
    output[idx] = delay * C # Convert to meter
    return output

def _sat_posvel(dset, time):
    file_key = "vlbi_orbit_sp3"
    rundate = dset.analysis["rundate"]
    days_before = (rundate - dset.time.datetime.min().date()).days
    days_after = (dset.time.datetime.max().date() - rundate).days
    orbit = apriori.get("basic_orbit", rundate=rundate,
                        file_key=file_key, days_before=days_before, days_after=days_after)

    idx = dset.near_field_obs
    num_sat_obs = np.sum(idx)

    satellites = np.unique(dset.source[idx])
    for sat in satellites:
        sat_pos = np.zeros((num_sat_obs, 3))
        sat_vel = np.zeros((num_sat_obs, 3))
        sat_idx = dset.source[idx] == sat
        sp3_sat_name = ngs_to_sp3[sat]
        sat_pos[sat_idx, :] = orbit[sp3_sat_name]["pos"](time[sat_idx])
        sat_vel[sat_idx, :] = orbit[sp3_sat_name]["vel"](time[sat_idx])
        sat_posvel = PosVel(np.concatenate((sat_pos, sat_vel), axis=1), system="trs", time=time)

    return sat_posvel

def _g2b_pos(r, t):
    """ Convert near Earth position expressed in GCRS to BRCS (TCG to TDB compatible).

    Based on equation 11.19 IERS2010 conventions and Deuv 2012 equation 4.

    Args:
    r:          Position in GCRS. Dimensions: (3, num_obs)
    t:          Time object for epoch for ephemerides (Length: num_obs)

    Returns:
    r_b:        Position in BCRS. Dimensions (3, num_obs)

    """
    eph = apriori.get("ephemerides", time=t)
    # The gravitational potential at the geocenter, neglecting the effects of the Earth’s mass.
    # At the picosecond level, only the solar potential is needed (IERS Conventions chapter 11 Table 11.1)
    U = constant.GM_sun / np.linalg.norm(eph.pos_gcrs("sun"), axis=1)[:, None]
    V_E = eph.vel_bcrs("earth")[:, None, :]
    R_E = eph.pos_bcrs("earth")
    #return R_E + r
    r_b = r * (1 - U / C ** 2 - L_C) - 0.5 * ((V_E @ r[:, :, None]) / C ** 2  @ V_E)[:, 0, :]
    return R_E + r_b

def _g2b_vel(v, t):
    """ Convert near Earth velocity expressed in GCRS to BCRS (TCG to TDB compatible)

    Based on Deuv 2012 equation 5

    Args:
    v:          Velocity in GCRS. Dimensions: (3, num_obs)
    t:          Time object for epoch of ephemerides (Length: num_obs)

    Returns:
    v_b:        Velocity in BCRS. Dimensions: (3, num_obs)
    """
    eph = apriori.get("ephemerides", time=t)
    U = constant.GM_sun / np.linalg.norm(eph.pos_gcrs("sun"), axis=1)[:, None]
    V_E = eph.vel_bcrs("earth")

    V_E_2 = (V_E[:, None, :] @ V_E[:, :, None])[:, :, 0]
    V_E_dot_v = (V_E[:, None, :] @ v[:, :, None])[:, :, 0]

    v_b =  (v * (1 - ((1 + GAMMA) * U)/ C ** 2 - V_E_2/(2 * C ** 2) - V_E_dot_v / C ** 2 ) 
            + V_E * (1 - 1/(2 * C ** 2) * V_E_dot_v))
    return v_b


def _deuv_relativistic_term(R_sat, R_site, T_sat, T_site, bodies):
    GM = _get_GM(bodies)

    eph_T_sat = apriori.get("ephemerides", time=T_sat)
    eph_T_site = apriori.get("ephemerides", time=T_site)
    
    R_sun_T_sat = eph_T_sat.pos_bcrs("sun")
    R_sun_T_site = eph_T_site.pos_bcrs("sun")
    
    R_sat_sun = R_sat - R_sun_T_sat # eq. 16, i = 0, alpha = sun
    R_site_sun = R_site - R_sun_T_site # eq. 16, i = 1 or 2, alpha = sun
    R_sat_site_sun = R_site_sun - R_sat_sun # eq. 17, alpha = sun

    norm_R_sat_sun = np.linalg.norm(R_sat_sun, axis=1)
    norm_R_site_sun = np.linalg.norm(R_site_sun, axis=1)
    norm_R_sat_site_sun = np.linalg.norm(R_sat_site_sun, axis=1)
    
    delay_sun = 0
    if "sun" in bodies:
        # The sun is treated individually in the following equations so remove it from the list of bodies
        bodies.remove("sun")

        # eq. 14 (first part)
        factor_sun = (1 + GAMMA) * GM["sun"] / C**3
        delay_sun = (factor_sun * C 
            * np.log((norm_R_sat_sun + norm_R_site_sun + norm_R_sat_site_sun + factor_sun) 
                     / (norm_R_sat_sun + norm_R_site_sun - norm_R_sat_site_sun + factor_sun)))

    delay_bodies = 0
    for body in bodies:
        
        R_body_T_sat = eph_T_sat.pos_bcrs(body)
        R_body_T_site = eph_T_site.pos_bcrs(body)
        
        R_sat_body = R_sat - R_body_T_sat # eq. 16, i = 0, alpha = body
        R_site_body = R_site  - R_body_T_site # eq. 16, i = 1 or 2, alpha = body
        R_sat_site_body = R_site_body - R_sat_body # eq. 17, alpha = body
        
        norm_R_sat_body = np.linalg.norm(R_sat_body, axis=1)
        norm_R_site_body = np.linalg.norm(R_site_body, axis=1)
        norm_R_sat_site_body = np.linalg.norm(R_sat_site_body, axis=1)
        
        # eq. 14 (last part)
        factor_body = (1 + GAMMA) * GM[body]/C ** 3
        delay_body = (factor_body 
            * np.log((norm_R_sat_body + norm_R_site_body + norm_R_sat_site_body)
                     / (norm_R_sat_body + norm_R_site_body - norm_R_sat_site_body)))
        delay_bodies += delay_body
        
        
    return delay_sun + delay_bodies # delay in TDB

@register_model
def deuv2012(dset):
    """
    Implementation based on 
        Deuv et al. 2012 : Spacecraft VLBI and Doppler tracking: algorithmns and implementations 
    """
    # This model is only applicable for near field observations
    idx = dset.near_field_obs
    num_sat_obs = np.sum(idx)
    if num_sat_obs == 0:
        # Skip this model if there are no near field observations
        return np.zeros(dset.num_obs)

    # This model is expressed in BCRS and should include all solar system bodies
    bodies = ["mercury", "venus", "earth", "moon", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto", "sun"]

    T_1 = dset.time.tdb[idx]
    r_1 = dset.site_pos_1.gcrs.pos.val[idx] # Site pos at t1 in GCRS
    R_1 = _g2b_pos(r_1, T_1)

    # Initial condition T0 = T1
    T_0 = T_1

    # Stop conditions
    limit = 0.1 * Unit.ps2s # TODO. When to break loop?
    max_iter = 5

    # Solve light-time equation for T_0 iteratively for signal path LT1 (from spacecraft to station 1)
    iter = 0
    while True:
        sat_posvel = _sat_posvel(dset, T_0) # Sat pos at t0 in ITRS

        r_0 = sat_posvel.gcrs.pos.val # Sat pos at t0 in GCRS
        v_0 = sat_posvel.gcrs.vel.val # Sat vel at t0 in GCRS
        R_0 = _g2b_pos(r_0, T_0) # Sat pos at T0 in BCRS
        V_0 = _g2b_vel(v_0, T_0) # Sat vel at T0 in BCRS

        RLT_01 = _deuv_relativistic_term(R_0, R_1, T_0, T_1, bodies) # Eq. 14
        R_01 = R_1 - R_0
        norm_R_01 = np.linalg.norm(R_01, axis=1)
        LT1 = norm_R_01 / C + RLT_01 # Eq. 13
        
        p_dot_01 = (R_01[:, None, :] / norm_R_01[:, None, None] @ V_0[:, :, None])[:, 0, 0] # Eq. 19
        delta_T_0 = (((T_1 - T_0).seconds - norm_R_01 / C - RLT_01) / 
                     (1 - p_dot_01 / C)) # Eq. 18
        iter += 1

        if any(np.abs(delta_T_0) < limit) or iter > max_iter:
            print(f"Stop after {iter-1} iterations")
            break

        #print(f"delta_T_0: {delta_T_0}")
        delta_T_0 = TimeDelta(delta_T_0, fmt="seconds", scale="tdb")
        T_0 = T_0 + delta_T_0

    # Initial condition for T2
    T_2 = T_1

    # Solve light time equation for T_2 iteratively for signal path LT2 (from spacecraft to station 2)
    iter = 0
    while True:
        # Assume that velocity of station 2 in trs is zero. Meaning trs pos is identical at T1 and T2
        station_2_trs = PosVel(dset.site_pos_2.val[idx], system="trs", time=T_2)

        r_2 = station_2_trs.gcrs.pos.val # Sat pos at t2 in GCRS
        v_2 = station_2_trs.gcrs.vel.val # Sat vel at t2 in GCRS
        R_2 = _g2b_pos(r_2, T_2) # Sat pos at T2 in BCRS
        V_2 = _g2b_vel(v_2, T_2) # Sat vel at T2 in BCRS

        RLT_02 = _deuv_relativistic_term(R_0, R_2, T_0, T_2, bodies)
        R_02 = R_2 - R_0
        norm_R_02 = np.linalg.norm(R_02, axis=1)
        LT2 = norm_R_02 / C + RLT_02

        p_dot_02 = (R_02[:, None, :] / norm_R_02[:, None, None] @ V_2[:, :, None])[:, 0, 0]
        delta_T_2 = (((T_2 - T_0).seconds - norm_R_02 / C - RLT_02) /
                     (-1 + p_dot_02 / C))
        iter += 1

        if any(np.abs(delta_T_2) < limit) or iter > max_iter:
            print(f"Stop after {iter-1} iterations")
            break

        #print(f"delta_T_2: {delta_T_2}")
        delta_T_2 = TimeDelta(delta_T_2, fmt="seconds", scale="tdb")
        T_2 = T_2 + delta_T_2

    delay_TDB = LT2 - LT1
    delay_TDB_2 = (T_2 - T_1).seconds
    #import IPython; IPython.embed()
    if any(np.abs(delay_TDB - delay_TDB_2)) > 1 * Unit.s2ps:
        print(f"(LT2 - LT1) - (T_2 - T_1) > 1 ps")

    # Convert to TT (Deuv 2012 eq. 20)
    eph = apriori.get("ephemerides", time=T_1)
    V_E = eph.vel_bcrs("earth")
    U = constant.GM_sun / np.linalg.norm(eph.pos_gcrs("sun"), axis=1) 
    b = (dset.site_pos_2.gcrs.pos.val[idx] - dset.site_pos_2.gcrs.pos.val[idx])[:, : , None]
    V_E_2 = (V_E[:, None, :] @ V_E[:, :, None])[:, 0, 0]
    V_E_dot_b = (V_E[:, None, :] @ b)[:, 0, 0]
    V_E_dot_v_2 = (V_E[:, None, :] @ v_2[:, :, None])[:, 0, 0]
    delay_TT = ((delay_TDB / (1 - L_C) * (1 - 1 / C ** 2 * (V_E_2 / 2 + U)) - V_E_dot_b / C ** 2) / 
                (1 + V_E_dot_v_2 / C ** 2))

    #import IPython; IPython.embed()
    _save_posvel_to_dset(dset, idx, f"{MODEL}.d_sat_posvel", sat_posvel, write_level="detail")
    _save_time_to_dset(dset, idx, f"{MODEL}.T_0", T_0, write_level="detail")
    _save_time_to_dset(dset, idx, f"{MODEL}.T_2", T_2, write_level="detail")
    _save_float_to_dset(dset, idx, f"{MODEL}.LT1", LT1, write_level="detail", unit="seconds")
    _save_float_to_dset(dset, idx, f"{MODEL}.LT2", LT2, write_level="detail", unit="seconds")
    _save_float_to_dset(dset, idx, f"{MODEL}.RLT_02", RLT_02, write_level="detail", unit="seconds")
    _save_float_to_dset(dset, idx, f"{MODEL}.RLT_01", RLT_01, write_level="detail", unit="seconds")

    output = np.zeros(dset.num_obs)
    output[idx] = delay_TT * C # Convert to meter
    return output

def _get_GM(bodies):
    GM = {}
    # Get GM for the celestial bodies
    for body in bodies:
        try:
            GM_name = "GM" if body == "earth" else f"GM_{body}"
            ephemerides = config.tech.ephemerides.str
            GM[body] = constant.get(GM_name, source=ephemerides)
        except KeyError:
            log.warn(
                f"The GM value of {body} is not defined for {ephemerides}. "
                f"Correction set to zero."
            )
            continue
    return GM

# Helper functions

def _save_float_to_dset(dset, idx, field, value, **kwargs):
    new_shape = tuple([dset.num_obs] + list(value.shape[1:]))
    full_value = np.full(new_shape, fill_value=np.nan)
    full_value[idx] = value
    dset.add_float(field, full_value, **kwargs)

def _save_posvel_to_dset(dset, idx, field, value, **kwargs):
    new_shape = tuple([dset.num_obs] + list(value.shape[1:]))
    full_value = np.full(new_shape, fill_value=np.nan)
    full_value[idx] = value.val
    dset.add_posvel(field, full_value, system=value.system, **kwargs)

def _save_time_to_dset(dset, idx, field, value, **kwargs):
    # Use datetime.min to indicate non-value
    t_min = Time(datetime.min, scale=value.scale, fmt="datetime")
    jd1 = np.full(dset.num_obs, fill_value=t_min.jd1)
    jd2 = np.full(dset.num_obs, fill_value=t_min.jd2)
    jd2[idx] = value.jd2
    jd1[idx] = value.jd1
    dset.add_time(field, val=jd1, val2=jd2, scale=value.scale, fmt="jd", **kwargs)

