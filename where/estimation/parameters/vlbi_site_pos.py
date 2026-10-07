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
# Standard library imports
from copy import deepcopy

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.constant import constant

# Where imports
from where.lib import config
from where.lib import rotation
from where.data.position import PosVel, PosVelDelta
from where.data.time import TimeDelta
from where.models.delay.vlbi_near_field import vlbi_near_field

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
        Tuple:    Array of partial derivatives, list of their names, and their unit
    """
    # Remove stations that should be fixed
    stations = np.asarray(dset.unique("station"))
    fix_stations = config.tech[PARAMETER].fix_stations.list
    fix_idx = np.in1d(stations, fix_stations)
    if fix_idx.any():
        stations = stations[np.logical_not(fix_idx)]

    idx = dset.near_field_obs
    nf_num_obs = np.sum(idx)

    column_names = [s + "_" + xyz for s in stations for xyz in "xyz"]
    partials = np.zeros((dset.num_obs, len(stations) * 3))
    partials_unit = "dimensionless"

    # Calculate partials for far field observations (typically quasars)
    all_partials = -dset.src_dir.unit_vector[:, None, :][~idx] @ rotation.trs2gcrs(dset.time)[~idx]

    # Calculate partials for near field observations (typically satellites)
    if nf_num_obs > 0:
        numerical_partials = config.tech[PARAMETER].numerical_partials.bool
        if numerical_partials:
            dtau_dx1, dtau_dx2 = _numerical_partials_near_field(dset)
        else:
            dtau_dx2 = _site_pos_2_near_field(dset)
            dtau_dx1 = _site_pos_1_near_field(dset)

    # Create partials matrix
    for i, station in enumerate(stations):
        filter_1 = dset.filter(station_1=station)
        partials[filter_1 & ~idx, i * 3 : i * 3 + 3] = all_partials[filter_1[~idx]][:, 0] * -1
        filter_2 = dset.filter(station_2=station)
        partials[filter_2 & ~idx, i * 3 : i * 3 + 3] = all_partials[filter_2[~idx]][:, 0]
        if nf_num_obs > 0:
            partials[filter_1 & idx, i * 3 : i * 3 + 3] = dtau_dx1[filter_1[idx]]
            partials[filter_2 & idx, i * 3 : i * 3 + 3] = dtau_dx2[filter_2[idx]]

    return partials, column_names, partials_unit

def  _site_pos_1_near_field(dset):
    """ Compute partial of near field observations with regards to station 1
    
    Args:
        data:     A Dataset containing model data
    
    Returns:
        np.array: Array with partial values. Dimensions (num_near_field_obs, 3)
    """
    idx = dset.near_field_obs
    num_obs = np.sum(idx)

    # Apriori values given at epoch t1 (time of arrival for signal at station 1)
    x1_t1 = dset.site_pos_1.gcrs.pos.val[idx] # station_1 at epoch t1
    x0_t1 = dset.sat_pos.gcrs.pos.val[idx] # satellite position at epoch t1
    v0_t1 = dset.sat_pos.gcrs.vel.val[idx] # satellite velocity at at epoch t1
    v2_t1 = dset.site_pos_2.gcrs.vel.val[idx] # station_2 velocity at epoch t1

    I = np.repeat(np.eye(3)[None, :,:], num_obs, axis=0)

    v0_t0_tilde = dset.vlbi_near_field.v0_t0_tilde[idx]
    v2_t2_tilde = dset.vlbi_near_field.v2_t2_tilde[idx]
    gamma0 = dset.vlbi_near_field.gamma0[idx, None, None]
    gamma2 = dset.vlbi_near_field.gamma2[idx, None, None]

    x1_x0 = x1_t1 - x0_t1
    ddelta_tilde_1_dx1 = x1_x0 / (C * np.linalg.norm(x1_x0, axis=1)[:, None]) # Eq 25 (skeens2024)
    dtau_tilde_dx1 = - ddelta_tilde_1_dx1 # Eq 32 (skeens2024) ***

    dx01_dx1 = v0_t0_tilde[:, : , None] @ ddelta_tilde_1_dx1[:, None, :] - I # Eq 24 (skeens2024)

    x01 = dset.vlbi_near_field.x01[idx]
    norm_x01 = np.linalg.norm(x01, axis=1)
    x01_dot_v0 = x01[:, None, :] @ v0_t1[:, :, None]  # Intermediate variable

    ddeltat0_dx1 = (1 / C ** 2 * dx01_dx1 @ v0_t1[:, :, None] 
        - gamma0 / C * (dx01_dx1 @ x01[:, :, None] + x01_dot_v0 * dx01_dx1 @ v0_t1[:, :, None] / C ** 2) 
        / np.sqrt(norm_x01 ** 2 + (x01_dot_v0[:, 0, 0]) ** 2 / C ** 2)[:, None, None]) # Eq 23 (skeens2024)

    dx02_dx1 = (v0_t0_tilde[:, :, None] @ ddelta_tilde_1_dx1[:, None, :] 
        + v2_t2_tilde[:, :, None] @ dtau_tilde_dx1[:, None, :]
        + (v0_t1 - v2_t1)[:, :, None] @ ddeltat0_dx1[:, :, 0][:, None, :]) # Eq 31 (skeens)

    x02 = dset.vlbi_near_field.x02[idx]
    norm_x02 = np.linalg.norm(x02, axis=1)
    x02_dot_v2 = x02[:, None, :] @ v2_t1[:, :, None]  # Intermediate variable

    ddeltat2_dx1 = (- 1 / C ** 2 * dx02_dx1 @ v2_t1[:, :, None]
        + gamma2 / C * (dx02_dx1 @ x02[:, :, None] + x02_dot_v2 * dx02_dx1 @ v2_t1[:, :, None] / C ** 2) 
        / np.sqrt(norm_x02 ** 2 + (x02_dot_v2[:, 0, 0]) ** 2 / C ** 2)[:, None, None]) # Eq 29 (skeens2024)

    # Do not convert to TCG-compatible values (skip (1 - L_G) from equation)
    dtau_dx1 = (ddeltat0_dx1 + ddeltat2_dx1) # Eq 22 (skeens)

    g2t = rotation.gcrs2trs(dset.time[idx])

    return (g2t @ dtau_dx1 * C)[:, :, 0] # Transform to trs and convert from seconds to meter

def _site_pos_2_near_field(dset):
    """ Compute partial of near field observations with regards to station 2
    
    Args:
        data:     A Dataset containing model data
    
    Returns:
        np.array: Array with partial values. Dimensions (num_near_field_obs, 3)
    """
    idx = dset.near_field_obs
    num_obs = np.sum(idx)

    # Apriori values given at epoch t1 (time of arrival for signal at station 1)
    x2_t1 = dset.site_pos_2.gcrs.pos[idx] # station_2 at epoch t1
    x0_t1 = dset.sat_pos.gcrs.pos[idx] # satellite position at epoch t1
    v2_t1 = dset.site_pos_2.gcrs.vel.val[idx] # station_2 velocity at epoch t1

    I = np.repeat(np.eye(3)[None, :,:], num_obs, axis=0)

    x2_x0 = x2_t1.val - x0_t1.val
    dtau_tilde_dx2 = x2_x0 / (C * np.linalg.norm(x2_x0, axis=1)[:, None]) # Eq 28 (skeens2024)

    v2_t2_tilde = dset.vlbi_near_field.v2_t2_tilde[idx]
    dx02_dx2 = v2_t2_tilde[:, :, None] @ dtau_tilde_dx2[:, None, :] - I # Eq 27 (skeens2024)

    x02 = dset.vlbi_near_field.x02[idx]
    norm_x02 = np.linalg.norm(x02, axis=1)
    x02_dot_v2 = (x02[:, None, :] @ v2_t1[:, :, None])[:, 0, 0]  # Intermediate variable

    gamma2 = dset.vlbi_near_field.gamma2[idx, None, None]

    ddeltat2_dx2 = (- 1 / C ** 2 * dx02_dx2 @ v2_t1[:, :, None]
        + (gamma2 / C * (dx02_dx2 @ x02[:, :, None] + x02_dot_v2[:, None, None] * dx02_dx2 @ v2_t1[:, :, None] / C ** 2))  
        / np.sqrt(norm_x02 ** 2 + (x02_dot_v2 / C) ** 2)[:, None, None]) # Eq 26 (skeens2024)

    # ddeltat0_dx2 = 0 according to eq 30 (skeens2024)

    # Do not convert to TCG-compatible values (skip (1 - L_G) from equation)
    dtau_dx2 = ddeltat2_dx2 # Eq 22 (skeens2024)

    g2t = rotation.gcrs2trs(dset.time[idx])
    return (g2t @ dtau_dx2  * C)[:, :, 0] # Transform to trs and convert from seconds to meter

def _numerical_partials_near_field(dset):
    """Compute numerical partial derivates for the near field delay.

    Uses the same near field model that was applied a priori.

    Args:
        dset                        A Dataset containing model data

    Returns:
        list [dtau_dx1, dtau_dx2]   Partial for station 1 and station 2. 
                                    Dimensions for each partial (num_near_field_obs x 3)
    """
  
    def _compute_numerical_partial(i):
        """Compute numerical partial for one coordinate component for one station

        Default field suffix determines which station

        Args:
            i (int)    0, 1 or 2 representing x, y, z

        Returns:
            dtau_dx_i
        """
        dx = np.zeros((dset.num_obs, 6))
        delta = 0.1
        dx[:, i] = delta
        dx = PosVelDelta(dx, system="trs", ref_pos=dset.site_pos)

        dset_1 = deepcopy(dset)
        dset_1.analysis = dset.analysis
        dset_1.default_field_suffix = dset.default_field_suffix
        dset_1.site_pos[:] -= dx
        d1 = vlbi_near_field(dset_1)

        dset_2 = deepcopy(dset)
        dset_2.analysis = dset.analysis
        dset_2.default_field_suffix = dset.default_field_suffix
        dset_2.site_pos[:] += dx
        d2 = vlbi_near_field(dset_2)

        return (d2 - d1) / (2 * delta)

    dtau_dx = []
    for _ in dset.for_each_suffix("station"):
        dtau_dx_i = np.zeros((dset.num_obs, 3))
        for i in range(3):
            dtau_dx_i[:, i] = _compute_numerical_partial(i)
        dtau_dx.append(dtau_dx_i)

    return dtau_dx

