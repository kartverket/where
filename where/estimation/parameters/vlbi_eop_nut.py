"""Calculate the partial derivatives of the celestial pole offset Earth Orientation Parameters.

Description:
------------

Calculate the partial derivatives of the :math:`X` and :math:`Y` Earth orientation parameters.

This is done according to equations (2.37) - (2.46) in Teke :cite:`teke2011`.




"""
# External library imports
import numpy as np

# Where imports
from where.lib import plugins
from where.lib import rotation


@plugins.register
def eop_nut(dset):
    """Calculate the partial derivative of the celestial pole offset Earth Orientation Parameters

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    column_names = ["x", "y"]
    partials = np.zeros((dset.num_obs, 2))
    src_dir = dset.src_dir.unit_vector[:, None, :]
    baseline = (dset.site_pos_2.itrs_pos - dset.site_pos_1.itrs_pos)[:, :, None]

    partials[:, 0] = -(src_dir @ rotation.dQ_dX(dset.time) @ rotation.R(dset.time) @ rotation.W(dset.time) @ baseline)[
        :, 0, 0
    ]
    partials[:, 1] = -(src_dir @ rotation.dQ_dY(dset.time) @ rotation.R(dset.time) @ rotation.W(dset.time) @ baseline)[
        :, 0, 0
    ]

    return partials, column_names, "meter per radian"
