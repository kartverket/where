"""Calculate the partial derivatives of the rate of the polar motion Earth Orientation Parameters.

Description:
------------

Calculate the partial derivatives of the rate of the :math:`x_p` and :math:`y_p` Earth orientation parameters.




$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""
# External library imports
import numpy as np

# Where imports
from where.ext import sofa_wrapper as sofa
from where.lib import plugins


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
    baseline = (dset.site_pos_2.itrs_pos - dset.site_pos_1.itrs_pos)[:, :, None]
    dt = (time.jd - time.mean.jd)[:, None, None]

    # x-pole
    partials[:, 0] = (src_dir @ sofa.Q(time) @ sofa.R(time) @ sofa.dW_dxp(time) @ baseline @ dt)[:, 0, 0]

    # y-pole
    partials[:, 1] = (src_dir @ sofa.Q(time) @ sofa.R(time) @ sofa.dW_dyp(time) @ baseline @ dt)[:, 0, 0]

    return partials, column_names, "meter * days / radian"
