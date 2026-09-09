"""Returns corrections in the geometric delay due to the propagation through the atmosphere.

Description:
------------

Calculate the geometric propagation delay using the Consensus model as described in the IERS Conventions
:cite:`iers2010`, section 11.1.




"""
# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from midgard.math.constant import constant
from where.lib import log


@plugins.register_ordered(1000)
def geometric_atmos(dset):
    """Returns the part of the geometric delay due to propagation through the atmosphere for each baseline

    This model depends on `troposphere_radio` already having run. Thus, the sort value is set to 1000 to make sure it
    runs last.

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array: Corrections in meters for each observation.

    """
    delay = np.zeros(dset.num_obs)
    idx = dset.near_field_obs
    if "troposphere_dT_1" in dset.fields:
        # Note that atm1 is already given in meter
        # The division by speed on light is a part of the model and not a unit conversion in this case
        atm1 = dset.troposphere_dT_1 / constant.c
    else:
        log.warn("Missing troposphere data. Atmospheric aberration correction set to zero")
        atm1 = np.zeros(dset.num_obs)

    # Far field model (from IERS 2010 Conventions)
    # Geometric delay due to the atmosphere in equation (11.11)
    baseline_gcrs_vel = (dset.site_pos_2.gcrs - dset.site_pos_1.gcrs).vel.val
    delay[~idx] = atm1[~idx] * (baseline_gcrs_vel[~idx][:, None, :] @ dset.src_dir.unit_vector[~idx][:, :, None])[:, 0, 0]

    # Near field model (from Hakan paper (unpublished))
    if np.sum(idx) > 0:
        k1 = dset.site_pos_1.gcrs.vector[idx]
        k1_hat = (k1 / np.linalg.norm(k1, axis=1)[:, None])[:, None, :]
        k2 = dset.site_pos_2.gcrs.vector[idx]
        k2_hat = k2 / np.linalg.norm(k2, axis=1)[:, None][:, None, :]
        v0 = dset.sat_pos.gcrs.vel.val[idx][:, :, None]
        v1 = dset.site_pos_1.gcrs.vel.val[idx][:, :, None]
        v2 = dset.site_pos_2.gcrs.vel.val[idx][:, :, None]
        delay[idx] = atm1[idx] * (k2_hat @ (v2 - v0) + k1_hat @ (v0 - v1))[:, 0, 0]

    # Since atm1 is already in meter we do not need to convert from seconds to meter
    return delay
