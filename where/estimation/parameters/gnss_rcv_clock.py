"""Calculate the partial derivatives of the GNSS receiver clock

Description:
------------

Calculate the partial derivatives of the GNSS receiver clock, e.g. described cite:`hofmann2008` p. 252.




"""
# External library imports
import numpy as np

# Where imports
from where.lib import constant
from where.lib import plugins


@plugins.register
def gnss_rcv_clock(dset):
    """Calculate the partial derivative of the GNSS receiver clock

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    station = dset.dataset_name
    partials = np.full((dset.num_obs, 1), constant.c)
    column_name = station + "_clock"

    return partials, column_name, "meter per second"
