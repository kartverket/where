"""Calculate the partial derivatives of the polar motion Earth Orientation Parameters.

Description:
------------

Calculate the partial derivatives of the :math:`x_p` and :math:`y_p` Earth orientation parameters.



"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import rotation


@plugins.register
def eop_pm(dset):
    """Calculate the partial derivative of the polar motion Earth Orientation Parameters

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    # Only far field observations is used to estimate this parameter
    idx = ~dset.near_field_obs

    column_names = ["xp", "yp"]
    partials = np.zeros((dset.num_obs, 2))

    time = dset.time[idx]
    src_dir = dset.src_dir.unit_vector[:, None, :][idx]
    baseline = (dset.site_pos_2.trs.pos[idx] - dset.site_pos_1.trs.pos[idx]).mat

    # x-pole
    partials[idx, 0] = -(
        src_dir @ rotation.Q(time) @ rotation.R(time) @ rotation.dW_dxp(time) @ baseline
    )[:, 0, 0]

    # y-pole
    partials[idx, 1] = -(
        src_dir @ rotation.Q(time) @ rotation.R(time) @ rotation.dW_dyp(time) @ baseline
    )[:, 0, 0]

    return partials, column_names, "meter per radian"
