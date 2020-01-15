"""Calculate the partial derivatives of the GNSS receiver clock

Description:
------------

Calculate the partial derivatives of the GNSS receiver clock, e.g. described cite:`hofmann2008` p. 252.

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.constant import constant


@plugins.register
def gnss_rcv_clock(dset):
    """Calculate the partial derivative of the GNSS receiver clock

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    partials = np.full((dset.num_obs, 1), constant.c)
    column_name = [""]

    return partials, column_name, "meter per second"
