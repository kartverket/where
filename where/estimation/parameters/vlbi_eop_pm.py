"""Calculate the partial derivatives of the polar motion Earth Orientation Parameters.

Description:
------------

Calculate the partial derivatives of the :math:`x_p` and :math:`y_p` Earth orientation parameters.

This is done according to equations (2.30) and (2.32) in Teke :cite:`teke2011`.




"""
# External library imports
import numpy as np

# Where imports
from where.ext import sofa_wrapper as sofa
from where.lib import plugins


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
    baseline = (dset.site_pos_2.itrs_pos - dset.site_pos_1.itrs_pos)[:, :, None]

    # x-pole
    partials[:, 0] = -(src_dir @ sofa.Q(dset.time) @ sofa.R(dset.time) @ sofa.dW_dxp(dset.time) @ baseline)[:, 0, 0]

    # y-pole
    partials[:, 1] = -(src_dir @ sofa.Q(dset.time) @ sofa.R(dset.time) @ sofa.dW_dyp(dset.time) @ baseline)[:, 0, 0]

    return partials, column_names, "meter per radian"
