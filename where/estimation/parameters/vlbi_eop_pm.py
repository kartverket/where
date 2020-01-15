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
    column_names = ["xp", "yp"]
    partials = np.zeros((dset.num_obs, 2))

    src_dir = dset.src_dir.unit_vector[:, None, :]
    baseline = (dset.site_pos_2.trs.pos - dset.site_pos_1.trs.pos).mat

    # x-pole
    partials[:, 0] = -(
        src_dir @ rotation.Q(dset.time) @ rotation.R(dset.time) @ rotation.dW_dxp(dset.time) @ baseline
    )[:, 0, 0]

    # y-pole
    partials[:, 1] = -(
        src_dir @ rotation.Q(dset.time) @ rotation.R(dset.time) @ rotation.dW_dyp(dset.time) @ baseline
    )[:, 0, 0]

    return partials, column_names, "meter per radian"
