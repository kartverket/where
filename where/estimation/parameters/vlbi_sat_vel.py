"""Calculate the partial derivatives of the satellite velocities

Description:
------------

Calculate the partial derivatives of the satellite velocities.

For near field targets the partial derivates is implemented based on the equations in cite:`skeens2024`


References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

.. [3] Skeens, Joe, Implementing a VLBI time delay model for Earth-orbiting satellites: partial derivaties and verification
       https://ntrs.nasa.gov/api/citations/20240007790/downloads/VTD%20Partials%20new%20update%20fmat.pdf


"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.constant import constant

# Where imports
from where.lib import config
from where.lib import rotation
from where.data.position import PosVel
from where.data.time import TimeDelta

# Name of parameter
PARAMETER = __name__.split(".")[-1]

# Constants for shorter equations
GAMMA = 1 # PPN parameter. Equal to 1 in general relativity
C = constant.c
L_G = constant.L_G
L_C = constant.L_C

@plugins.register
def site_pos(dset):
    """Calculate the partial derivative of the satellite position for each satellite

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple:    Array of partial derivatives, list of their names, and their unit
    """
    # Remove stations that should be fixed
    # stations = np.asarray(dset.unique("station"))
    # fix_stations = config.tech[PARAMETER].fix_stations.list
    # fix_idx = np.in1d(stations, fix_stations)
    # if fix_idx.any():
    #     stations = stations[np.logical_not(fix_idx)]

    idx = dset.near_field_obs
    nf_num_obs = np.sum(idx)
    
    satellites = np.unique(dset.source[dset.near_field_obs])

    # Calculate partials for near field observations (typically satellites)
    if nf_num_obs > 0:
        dtau_dv0 = _sat_vel(dset)
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots(3, sharex=True); label = "xyz";
    # for i in range(3):
    #     ax[i].scatter(dset.time.mjd[idx], dtau_dx1[:, i, 0], marker=".", alpha=0.5, label="dtau_dx1")
    #     ax[i].scatter(dset.time.mjd[idx], dtau_dx2[:, i, 0], marker=".", alpha=0.5, label="dtau_dx2")
    #     ax[i].scatter(dset.time.mjd[~idx], all_partials[:, 0, i], marker=".", alpha=0.5, label="dtau_dx")
    #     ax[i].set_ylabel(label[i])
    # plt.xlabel("mjd"); ax[0].legend(); plt.show()

    partials = np.zeros((dset.num_obs, len(satellites) * 3))
    for i, sat in enumerate(satellites):
        filter = dset.filter(source=sat)
        #partials[filter_1 & ~idx, i * 3 : i * 3 + 3] = all_partials[filter_1[~idx]][:, 0] * -1
        #filter_2 = dset.filter(source=sat)
        #partials[filter_2 & ~idx, i * 3 : i * 3 + 3] = all_partials[filter_2[~idx]][:, 0]
        if nf_num_obs > 0:
            partials[filter & idx, i * 3 : i * 3 + 3] = dtau_dv0[filter[idx]][:, :, 0]
            #partials[filter_2 & idx, i * 3 : i * 3 + 3] = dtau_dx2[filter_2[idx]][:, :, 0]

    column_names = [s + "_" + xyz for s in satellites for xyz in ["vx", "vy", "vz"]]

    return partials, column_names, "seconds"

def  _sat_vel(dset):
    """ Compute partial of near field observations with regards to the satellite
    
    Args:
        data:     A Dataset containing model data
    
    Returns:
        np.array: Array with partial values. Dimensions (num_near_field_obs, 3, 1)
    """
    idx = dset.near_field_obs
    num_obs = np.sum(idx)

    # Apriori values given at epoch t1 (time of arrival for signal at station 1)
    x1_t1 = dset.site_pos_1.gcrs.pos[idx] # station_1 at epoch t1
    x0_t1 = dset.sat_pos.gcrs.pos[idx] # satellite position at epoch t1
    v0_t1 = dset.sat_pos.gcrs.vel.val[idx] # satellite velocity at at epoch t1
    x2_t1 = dset.site_pos_2.gcrs.pos[idx] # station_2 at epoch t1
    v2_t1 = dset.site_pos_2.gcrs.vel.val[idx] # station_2 velocity at epoch t1

    I = np.repeat(np.eye(3)[None, :,:], num_obs, axis=0)

    v0_t0_tilde = dset.vlbi_near_field.v0_t0_tilde[idx]
    v2_t2_tilde = dset.vlbi_near_field.v2_t2_tilde[idx]
    gamma0 = dset.vlbi_near_field.gamma0[idx, None, None]
    gamma2 = dset.vlbi_near_field.gamma2[idx, None, None]

    # Since dx01_dv0 = 0 the equations below have been simplified
    #dx01_dv0 = 0 # Eq 59 (skeens2024)

    v0_2 = v0_t1[:, None, :] @ v0_t1[:, :, None]
    denominator = C ** 2 * (1 - v0_2 / C ** 2) ** (3/2)
    dgamma0_dv0 = v0_t1[:, :, None] / denominator # Eq 60 (skeens2024)
    
    x01 = dset.vlbi_near_field.x01[idx]
    norm_x01 = np.linalg.norm(x01, axis=1)
    x01_dot_v0 = x01[:, None, :] @ v0_t1[:, :, None]  # Intermediate variable
    sqrt_x01_v0 = np.sqrt(norm_x01 ** 2 + (x01_dot_v0[:, 0, 0]) ** 2 / C ** 2)[:, None, None] # Intermediate variable
    
    # TODO: test dimensions
    
    ddeltat0_dv0 = (1 / C ** 2 * x01[:, :, None]
        - 1 / C * dgamma0_dv0 * sqrt_x01_v0 
        - gamma0 / C * (x01_dot_v0 * x01[:, :, None] / C ** 2) 
        / sqrt_x01_v0) # Eq 58 simplified (skeens2024)

    delta_t0 = dset.vlbi_near_field.delta_t0[idx]
    
    dx02_dv0 = (delta_t0[:, None, None] * I
        + (v0_t1 - v2_t1)[:, :, None] @ ddeltat0_dv0[:, :, 0][:, None, :]) # Eq 62 (skeens)

    x02 = dset.vlbi_near_field.x02[idx]
    norm_x02 = np.linalg.norm(x02, axis=1)
    x02_dot_v2 = x02[:, None, :] @ v2_t1[:, :, None]  # Intermediate variable

    ddeltat2_dv0 = (- 1 / C ** 2 * dx02_dv0 @ v2_t1[:, :, None]
        + gamma2 / C * (dx02_dv0 @ x02[:, :, None] + x02_dot_v2 * dx02_dv0 @ v2_t1[:, :, None] / C ** 2) 
        / np.sqrt(norm_x02 ** 2 + (x02_dot_v2)[:, 0, 0] ** 2 / C ** 2)[:, None, None]) # Eq 61 (skeens2024)

    dtau_dv0 = (ddeltat0_dv0 + ddeltat2_dv0) * (1 - L_G) # Eq 50 (skeens)

    return dtau_dv0 * C # Convert from seconds to meter
