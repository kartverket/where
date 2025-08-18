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
from where.lib import log

GAMMA = 1 # PPN parameter. Equal to 1 in general relativity

@plugins.register
def vlbi_near_field_delay(dset):
    r"""Calculate the theoretical delay dependent on the baseline
    -------------------------------------------------------

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array: Projected baseline in meters for each observation.

    """
    orbit = apriori.get("basic_orbit", rundate=dset.analysis["rundate"], days_before=0, days_after=1)
    eph = apriori.get("ephemerides", time=dset.time)
    bodies = [
        "mercury",
        "venus",
        "earth",
        "moon",
        "mars",
        "jupiter",
        "saturn",
        "uranus",
        "neptune",
        "pluto",
        "sun",
    ]
    GM = {}
    # Get GM for the celestial bodies
    for body in bodies:
        try:
            GM_name = "GM" if body == "earth" else f"GM_{body}"
            GM[body] = constant.get(GM_name, source=eph.ephemerides)
        except KeyError:
            log.warn(
                f"The GM value of {body} is not defined for {eph.ephemerides}. "
                f"Correction set to zero."
            )
            continue
    # The sun is treated individually in the following equations so remove it from the list of bodies
    bodies.remove("sun")
    # TODO
    # idx_sat = True when observation is to a satellite

    delay = np.zeros(dset.num_obs)
    
    # Apriori values given at epoch t1 (time of arrival for signal at station 1)
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
    
    t0_tilde = t1 - delta1 # approximation to t0 (time of emission of signal from satellite)
    tau_tilde = delta2 - delta1 # eq. 7
    t2_tilde = t1 + tau_tilde # approximation to t2 (time of arrival for signal at station 2)
 
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
    

    # Compute t_g01: Relativistic effects on delay from satellite to station 1
    # Based on Deuv, et al (2012) eq. 14, 16, 17
    # Equations are in BCRS. Ephemerides use TDB.
    # Ref. IERS 2010 conventions eq. 11.6 used for grav delay in consensus model
    # Transformation between BCRS and GRCS is approximated by X_bcrs = X_earth + X_gcrs
    R0_T0 = eph.pos_bcrs("earth", time=t0_tilde) + x0_t0_tilde # satellite position at t0 in BCRS
    RS_T0 = eph.pos_bcrs("sun", time=t0_tilde) # sun pos at t0 in BCRS
    R1_T1 = eph.pos_bcrs("earth") + dset.site_pos_1.gcrs.pos.val # station_1 pos at t1 in BCRS
    RS_T1 = eph.pos_bcrs("sun") # sun pos at t1 in BCRS
    R0_S = R0_T0 - RS_T0 # eq. 16, i=0, alpha = S
    R1_S = R1_T1 - RS_T1 # eq. 16, i=1, alpha = S
    R01_S = R1_S - R0_S # eq. 17, alpha = S

    norm_R0_S = np.linalg.norm(R0_S, axis=1)
    norm_R1_S = np.linalg.norm(R1_S, axis=1)
    norm_R01_S = np.linalg.norm(R01_S, axis=1)

    # eq. 14 (first part)
    sun_factor = (1 + GAMMA) * GM["sun"]/constant.c ** 2
    delay_sun = sun_factor/constant.c * \
        np.log((norm_R0_S + norm_R1_S + norm_R01_S + sun_factor)/(norm_R0_S + norm_R1_S - norm_R01_S + sun_factor)) 

    _save_detail_to_dataset(dset, "vlbi_nf_grav_sun_1", delay_sun * constant.c, dset.add_float, unit="meter")

    delay_bodies = 0
    for body in bodies:
        RB_T0 = eph.pos_bcrs(body, time=t0_tilde) # body pos at t0 in BCRS
        RB_T1 = eph.pos_bcrs(body) # body pos at t1 in BCRS
        R0_B = R0_T0 - RB_T0 # eq. 16, i=0, alpha = B
        R1_B = R1_T1 - RB_T1 # eq. 16, i=1, alpha = B
        R01_B = R1_B - R0_B # eq. 17, aplha = B

        norm_R0_B = np.linalg.norm(R0_B, axis=1)
        norm_R1_B = np.linalg.norm(R1_B, axis=1)
        norm_R01_B = np.linalg.norm(R01_B, axis=1)

        # eq. 14 (last part)
        factor_body = (1 + GAMMA) * GM[body]/constant.c ** 3
        delay_body = factor_body * \
            np.log((norm_R0_B + norm_R1_B + norm_R01_B)/(norm_R0_B + norm_R1_B - norm_R01_B))
        delay_bodies += delay_body
        
        _save_detail_to_dataset(dset, f"vlbi_nf_grav_{body}_1", delay_body * constant.c, dset.add_float, unit="meter")

    #t_g01_TDB = np.zeros(dset.num_obs) 
    t_g01_TDB = delay_sun + delay_bodies # eq. 14 in deuv2012
    # According to Kaplan 2005: "TDB advance, on average, at the same rate as TT". 
    # -> Assume delay in TDB is the same as the delay in TT for this purpose
    # Convert from TT to TCG since the Jaron, et. al (2017) equations work with this
    t_g01 = t_g01_TDB / (1 - constant.L_G)
    
    # Save TT(=TDB) value to dset
    _save_detail_to_dataset(dset, "vlbi_nf_grav_1", t_g01_TDB * constant.c, dset.add_float, unit="meter")

    #import IPython; IPython.embed()
    
    # eq. 14 in jaron2017
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
    
    # Compute t_g02: Relativistic effects on delay from satellite to station 2
    # Based on Deuv, et al (2012) eq. 14, 16, 17
    # Equations are in BCRS. Ephemerides use TDB.  

    # Ref. IERS 2010 conventions eq. 11.6 used for grav delay in consensus model
    # Transformation between BCRS and GRCS is approximated by X_bcrs = X_earth + X_gcrs
    R2_T2 = eph.pos_bcrs("earth", time=t2_tilde) + x2_t2_tilde # station_2 pos at t2 in BCRS
    RS_T2 = eph.pos_bcrs("sun", time=t2_tilde) # sun pos at t2 in BCRS
    R2_S = R2_T2 - RS_T2 # eq. 16, i=2, alpha = S
    R02_S = R2_S - R0_S # eq. 17, alpha = S

    norm_R2_S = np.linalg.norm(R2_S, axis=1)
    norm_R02_S = np.linalg.norm(R02_S, axis=1)

    # eq. 14 (first part)
    sun_factor = (1 + GAMMA) * GM["sun"]/constant.c ** 2
    delay_sun = sun_factor/constant.c * \
        np.log((norm_R0_S + norm_R2_S + norm_R02_S + sun_factor)/(norm_R0_S + norm_R2_S - norm_R02_S + sun_factor)) 

    _save_detail_to_dataset(dset, "vlbi_nf_grav_sun_2", delay_sun * constant.c, dset.add_float, unit="meter")

    delay_bodies = 0
    for body in bodies:
        RB_T0 = eph.pos_bcrs(body) # body pos at t0 in BCRS
        RB_T2 = eph.pos_bcrs(body, time=t2_tilde) # body pos at t2 in BCRS
        R0_B = R0_T0 - RB_T0 # eq. 16, i=0, alpha = B
        R2_B = R2_T2 - RB_T2 # eq. 16, i=2, alpha = B
        R02_B = R2_B - R0_B # eq. 17, aplha = B

        norm_R0_B = np.linalg.norm(R0_B, axis=1)
        norm_R2_B = np.linalg.norm(R2_B, axis=1)
        norm_R02_B = np.linalg.norm(R02_B, axis=1)

        # eq. 14 (last part)
        factor_body = (1 + GAMMA) * GM[body]/constant.c ** 3
        delay_body = factor_body * \
            np.log((norm_R0_B + norm_R2_B + norm_R02_B)/(norm_R0_B + norm_R2_B - norm_R02_B))
        delay_bodies += delay_body 
        
        _save_detail_to_dataset(dset, f"vlbi_nf_grav_{body}_2", delay_body * constant.c, dset.add_float, unit="meter")

    #t_g02_TDB = np.zeros(dset.num_obs) 
    t_g02_TDB = delay_sun + delay_bodies # eq. 14 in deuv2012
    # According to Kaplan 2005: "TDB advance, on average, at the same rate as TT". 
    # -> Assume delay in TDB is the same as the delay in TT for this purpose
    # Convert from TT to TCG since the Jaron, et. al (2017) equations work with this
    t_g02 = t_g02_TDB / (1 - constant.L_G)
    
    # Save TT(=TDB) value to dset  
    _save_detail_to_dataset(dset, "vlbi_nf_grav_2", t_g02_TDB * constant.c, dset.add_float, unit="meter")
    
    # eq. 17 in jaron2017
    x_dot_v_2 = (x02[:, None, :] @ v2_t1[:, :, None])[:, 0, 0] / constant.c ** 2 # Intermediate variable
    x02_dot_x02 = (x02[:, None, :] @ x02[:, :, None])[:, 0, 0] # Intermediate variable
    # Time of reception of the signal relative to t1
    delta_t2 = - gamma2_2 * (x_dot_v_2 - t_g02) + \
        np.sqrt(gamma2_2 ** 2 * (x_dot_v_2 - t_g02) ** 2 + gamma2_2 * (x02_dot_x02 / constant.c ** 2 - t_g02 ** 2))
    
    # Convert from TCG to TT
    delay = (delta_t2 + delta_t0) * (1 - constant.L_G) # eq. 10 


    
    ## For debugging. See if satellite is above horizon for both stations
    s1 = dset.site_pos_1.copy()
    s1.other = dset.sat_pos
    s2 = dset.site_pos_2.copy()
    s2.other = dset.sat_pos
    sat_visible = (s1.elevation > 0) & (s2.elevation > 0)
    
    _save_detail_to_dataset(dset, "sat_visible", sat_visible, dset.add_bool)

    import IPython; IPython.embed()

    return delay * constant.c # Convert to meter

def _save_detail_to_dataset(dset, field, value, func, **kwargs):
    if field in dset.fields:
        dset[field][:] = value
    else:
        func(field, value, write_level="detail", **kwargs)
