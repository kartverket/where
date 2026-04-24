"""Calculate the partial derivatives of the site positions

Description:
------------

Calculate the partial derivatives of the site positions.

For far field targets the partial derivates is implemented based on the equations in cite:`teke2011`

For near field targets the partial derivates is implemented based on the equations in cite:`skeens2024`


References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

.. [2] Teke, Kamil, Sub-daily parameter estimation in VLBI data analysis.
       https://geo.tuwien.ac.at/fileadmin/editors/GM/GM87_teke.pdf

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
    """Calculate the partial derivative of the site position for each station

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, list of their names, and their unit
    """
    # Remove stations that should be fixed
    stations = np.asarray(dset.unique("station"))
    fix_stations = config.tech[PARAMETER].fix_stations.list
    fix_idx = np.in1d(stations, fix_stations)
    if fix_idx.any():
        stations = stations[np.logical_not(fix_idx)]

    #import IPython; IPython.embed()

    # Calculate partials for far field observations (typically quasars)
    all_partials = -dset.src_dir.unit_vector[:, None, :] @ rotation.trs2gcrs(dset.time)
    # Calculate partials for near field observations (typically satellites)
    # TODO
    # TODO if vlbi_near_field not in dset and dset.near_field_obs > 0 -> Error
    dtau_dx2 = _site_pos_2_near_field(dset)
    dtau_dx1 = _site_pos_1_near_field(dset)
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, sharex=True); label = "xyz";
    for i in range(3):
        ax[i].scatter(dset.time.mjd, dtau_dx1[:, i, 0], label="dtau_dx1")
        ax[i].scatter(dset.time.mjd, dtau_dx2[:, i, 0], label="dtau_dx2")
        ax[i].scatter(dset.time.mjd, all_partials[:, 0, i], label="dtau_dx")
        ax[i].set_ylabel(label[i])
    plt.xlabel("mjd"); ax[0].legend(); plt.show()

    import IPython; IPython.embed()
    all_partials[dset.near_field_obs] = 0
    partials = np.zeros((dset.num_obs, len(stations) * 3))
    for idx, station in enumerate(stations):
        filter_1 = dset.filter(station_1=station)
        partials[filter_1, idx * 3 : idx * 3 + 3] = all_partials[filter_1][:, 0] * -1
        filter_2 = dset.filter(station_2=station)
        partials[filter_2, idx * 3 : idx * 3 + 3] = all_partials[filter_2][:, 0]

    column_names = [s + "_" + xyz for s in stations for xyz in "xyz"]

    return partials, column_names, "dimensionless"

def  _site_pos_1_near_field(dset):
    partials = 0
    # Apriori values given at epoch t1 (time of arrival for signal at station 1)
    x1_t1 = dset.site_pos_1.gcrs.pos # station_1 at epoch t1
    x0_t1 = dset.sat_pos.gcrs.pos # satellite position at epoch t1
    v0_t1 = dset.sat_pos.gcrs.vel.val # satellite velocity at at epoch t1
    v2_t1 = dset.site_pos_2.gcrs.vel.val # station_2 velocity at epoch t1
    
    I = np.repeat(np.eye(3)[None, :,:], dset.num_obs, axis=0)
    
    v0_t0_tilde = dset.vlbi_near_field.v0_t0_tilde
    v2_t2_tilde = dset.vlbi_near_field.v2_t2_tilde
    gamma0 = dset.vlbi_near_field.gamma0[:, None, None]
    gamma2 = dset.vlbi_near_field.gamma2[:, None, None]
    
    x1_x0 = x1_t1.val - x0_t1.val
    ddelta_tilde_1_dx1 = (x1_x0) / (C * np.linalg.norm(x1_x0)) # Eq 25 (skeens2024)
    dtau_tilde_dx1 = - ddelta_tilde_1_dx1 # Eq 32 (skeens2024)
    
    dx01_dx1 = v0_t0_tilde[:, : , None] @ ddelta_tilde_1_dx1[:, None, :] - I # Eq 24 (skeens2024)
    
    x01 = dset.vlbi_near_field.x01
    norm_x01 = np.linalg.norm(x01)
    x01_dot_v0 = x01[:, None, :] @ v0_t1[:, :, None]  # Intermediate variable
    
    ddeltat0_dx1 = (1 / C ** 2 * dx01_dx1 @ v0_t1[:, :, None] 
        - gamma0 / C * (dx01_dx1 @ x01[:, :, None] + x01_dot_v0 * dx01_dx1 @ v0_t1[:, :, None] / C ** 2) 
        / np.sqrt(norm_x01 ** 2 + (x01_dot_v0[:, 0, 0]) ** 2 / C ** 2)[:, None, None]) # Eq 23 (skeens2024)
    
    dx02_dx1 = (v0_t0_tilde[:, :, None] @ ddelta_tilde_1_dx1[:, None, :] 
        + v2_t2_tilde[:, :, None] @ dtau_tilde_dx1[:, None, :]
        + (v0_t1 - v2_t1)[:, :, None] @ ddeltat0_dx1[:, :, 0][:, None, :]) # Eq 31 (skeens)
    
    x02 = dset.vlbi_near_field.x02
    norm_x02 = np.linalg.norm(x02)
    x02_dot_v2 = x02[:, None, :] @ v2_t1[:, :, None]  # Intermediate variable
    
    ddeltat2_dx1 = (- 1 / C ** 2 * dx02_dx1 @ v2_t1[:, :, None]
        + gamma2 / C * (dx02_dx1 @ x02[:, :, None] + x02_dot_v2 * dx02_dx1 @ v2_t1[:, :, None] / C ** 2) 
        / np.sqrt(norm_x02 ** 2 + (x02_dot_v2)[:, 0, 0] ** 2 / C ** 2)[:, None, None]) # Eq 29 (skeens2024)
    
    dtau_dx1 = (ddeltat0_dx1 + ddeltat2_dx1) * (1 - L_G) # Eq 22 (skeens)
    return dtau_dx1 * C # Convert from seconds to meter

def _site_pos_2_near_field(dset):
    
    partials = 0
    #filter_2 = dset.filter(station_2=station)
    
    # Apriori values given at epoch t1 (time of arrival for signal at station 1)
    x2_t1 = dset.site_pos_2.gcrs.pos # station_2 at epoch t1
    x0_t1 = dset.sat_pos.gcrs.pos # satellite position at epoch t1
    v2_t1 = dset.site_pos_2.gcrs.vel.val # station_2 velocity at epoch t1
    
    I = np.repeat(np.eye(3)[None, :,:], dset.num_obs, axis=0)
    
    x2_x0 = x2_t1.val - x0_t1.val
    dtau_tilde_dx2 = x2_x0 / (C * np.linalg.norm(x2_x0)) # Eq 28 (skeens2024)
    
    v2_t2_tilde = dset.vlbi_near_field.v2_t2_tilde
    dx02_dx2 = v2_t2_tilde[:, :, None] @ dtau_tilde_dx2[:, None, :] - I # Eq 27 (skeens2024)
    
    x02 = dset.vlbi_near_field.x02
    norm_x02 = np.linalg.norm(x02)
    x02_dot_v2 = (x02[:, None, :] @ v2_t1[:, :, None])[:, 0, 0]  # Intermediate variable
    
    gamma2 = dset.vlbi_near_field.gamma2[:, None, None]
    #import IPython; IPython.embed()
    ddeltat2_dx2 = (- 1 / C ** 2 * dx02_dx2 @ v2_t1[:, :, None]
        + (gamma2 / C * (dx02_dx2 @ x02[:, :, None] + x02_dot_v2[:, None, None] * dx02_dx2 @ v2_t1[:, :, None] / C ** 2))  
        / np.sqrt(norm_x02 ** 2 + (x02_dot_v2 / C) ** 2)[:, None, None]) # Eq 26 (skeens2024)
    
    # ddeltat0_dx2 = 0 according to eq 30 (skeens)
    dtau_dx2 = ddeltat2_dx2 * (1 - L_G) # Eq 22 (skeens2024)
    return dtau_dx2 * C # Convert from seconds to meter
