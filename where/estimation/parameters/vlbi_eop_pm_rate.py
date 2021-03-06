"""Calculate the partial derivatives of the rate of the polar motion Earth Orientation Parameters.

Description:
------------

Calculate the partial derivatives of the rate of the :math:`x_p` and :math:`y_p` Earth orientation parameters.





"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import rotation


@plugins.register
def eop_pm_rate(dset):
    """Calculate the partial derivative of the rate of the polar motion Earth Orientation Parameters

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    column_names = ["dxp", "dyp"]
    partials = np.zeros((dset.num_obs, 2))

    time = dset.time
    src_dir = dset.src_dir.unit_vector[:, None, :]
    baseline = (dset.site_pos_2.trs.pos - dset.site_pos_1.trs.pos).mat
    dt = (time.jd - time.mean.jd)[:, None, None]

    # x-pole
    partials[:, 0] = -(src_dir @ rotation.Q(time) @ rotation.R(time) @ rotation.dW_dxp(time) @ baseline @ dt)[:, 0, 0]

    # y-pole
    partials[:, 1] = -(src_dir @ rotation.Q(time) @ rotation.R(time) @ rotation.dW_dyp(time) @ baseline @ dt)[:, 0, 0]

    return partials, column_names, "meter * days / radian"
