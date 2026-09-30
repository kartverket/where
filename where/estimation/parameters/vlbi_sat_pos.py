"""Calculate the partial derivatives of the satellite positions

Description:
------------

Calculate the partial derivatives of the satellite positions.

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
def sat_pos(dset):
    """Calculate the partial derivative of the satellite position for each satellite

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple:    Array of partial derivatives, list of their names, and their unit
    """
    idx = dset.near_field_obs
    satellites = np.unique(dset.source[dset.near_field_obs])

    column_names = [s + "_" + xyz for s in satellites for xyz in "xyz"]
    partials = partials = np.zeros((dset.num_obs, len(satellites) * 3))
    partials_unit = "dimensionless"

    if np.sum(idx) == 0:
        return partials, column_names, partials_unit

    # Calculate partials for near field observations (typically satellites)
    dtau_dx0 = _sat_pos(dset)

    for i, sat in enumerate(satellites):
        filter = dset.filter(source=sat)
        partials[filter & idx, i * 3 : i * 3 + 3] = dtau_dx0[filter[idx]][:, :, 0]

    return partials, column_names, partials_unit

def  _sat_pos(dset):
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

    x1_x0 = x1_t1.val - x0_t1.val
    ddelta_tilde_1_dx0 = (x1_x0) / (C * np.linalg.norm(x1_x0, axis=1)[:, None]) # Eq 54 (skeens2024)
    
    x2_x0 = x2_t1.val - x0_t1.val
    ddelta_tilde_2_dx0 = -(x2_x0) / (C * np.linalg.norm(x2_x0, axis=1)[:, None]) # Eq 57 (skeens2024)
    dtau_tilde_dx0 = ddelta_tilde_2_dx0 - ddelta_tilde_1_dx0 # Eq 57 (skeens2024)

    dx01_dx0 = v0_t0_tilde[:, : , None] @ ddelta_tilde_1_dx0[:, None, :] + I # Eq 53 (skeens2024)

    x01 = dset.vlbi_near_field.x01[idx]
    norm_x01 = np.linalg.norm(x01, axis=1)
    x01_dot_v0 = x01[:, None, :] @ v0_t1[:, :, None]  # Intermediate variable

    ddeltat0_dx0 = (1 / C ** 2 * dx01_dx0 @ v0_t1[:, :, None] 
        - gamma0 / C * (dx01_dx0 @ x01[:, :, None] + x01_dot_v0 * dx01_dx0 @ v0_t1[:, :, None] / C ** 2) 
        / np.sqrt(norm_x01 ** 2 + (x01_dot_v0[:, 0, 0]) ** 2 / C ** 2)[:, None, None]) # Eq 52 (skeens2024)

    dx02_dx0 = (v0_t0_tilde[:, :, None] @ ddelta_tilde_1_dx0[:, None, :]
        + I 
        + v2_t2_tilde[:, :, None] @ dtau_tilde_dx0[:, None, :]
        + (v0_t1 - v2_t1)[:, :, None] @ ddeltat0_dx0[:, :, 0][:, None, :]) # Eq 56 (skeens)

    x02 = dset.vlbi_near_field.x02[idx]
    norm_x02 = np.linalg.norm(x02, axis=1)
    x02_dot_v2 = x02[:, None, :] @ v2_t1[:, :, None]  # Intermediate variable

    ddeltat2_dx0 = (- 1 / C ** 2 * dx02_dx0 @ v2_t1[:, :, None]
        + gamma2 / C * (dx02_dx0 @ x02[:, :, None] + x02_dot_v2 * dx02_dx0 @ v2_t1[:, :, None] / C ** 2) 
        / np.sqrt(norm_x02 ** 2 + (x02_dot_v2)[:, 0, 0] ** 2 / C ** 2)[:, None, None]) # Eq 55 (skeens2024)

    dtau_dx0 = (ddeltat0_dx0 + ddeltat2_dx0) * (1 - L_G) # Eq 48 (skeens)

    return dtau_dx0 * C # Convert from seconds to meter
