"""Calculate the partial derivatives of the celestial pole offset Earth Orientation Parameters.

Description:
------------

Calculate the partial derivatives of the :math:`X` and :math:`Y` Earth orientation parameters.


"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import rotation


@plugins.register
def eop_nut(dset):
    """Calculate the partial derivative of the celestial pole offset Earth Orientation Parameters

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    # Only far field observations is used to estimate this parameter
    idx = ~dset.near_field_obs

    column_names = ["x", "y"]
    partials = np.zeros((dset.num_obs, 2))

    time = dset.time[idx]   
    src_dir = dset.src_dir.unit_vector[:, None, :][idx]
    baseline = (dset.site_pos_2.trs.pos[idx] - dset.site_pos_1.trs.pos[idx]).mat

    partials[idx, 0] = -(src_dir @ rotation.dQ_dX(time) @ rotation.R(time) @ rotation.W(time) @ baseline)[
        :, 0, 0
    ]
    partials[idx, 1] = -(src_dir @ rotation.dQ_dY(time) @ rotation.R(time) @ rotation.W(time) @ baseline)[
        :, 0, 0
    ]

    return partials, column_names, "meter per radian"
