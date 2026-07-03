"""Calculate the partial derivatives of the vlbi scale

Description:
------------

Calculate the partial derivatives of the vlbi scale.

Implementation is based on cite:`titov2018` equation 10.

"""

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import rotation

# Name of parameter
PARAMETER = __name__.split(".")[-1]


@plugins.register
def scale(dset):
    """Calculate the partial derivative of the vlbi scale

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, list of their names, and their unit
    """
    # Only far field observations are used to estimate this parameter (for now)
    idx = ~dset.near_field_obs
    partials = np.zeros((dset.num_obs, 1))

    src_dir = dset.src_dir[idx].unit_vector[:, None, :]
    baseline = (dset.site_pos_2[idx].trs.pos - dset.site_pos_1[idx].trs.pos).mat
    partials[idx, :] = -(src_dir @ rotation.trs2gcrs(dset.time[idx]) @ baseline)[:, :, 0]

    return partials, ["scale"], "meter"
