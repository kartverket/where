"""Returns delay correction due to the cables

Description:
------------

The cable calibration is already calculated by the correlators and is provided on the NGS file.

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def cable_calibration(dset):
    """Calculate total delay due to cable calibration

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing corrections in meters for each observation
    """
    return -(np.nan_to_num(dset.cable_delay_2) - np.nan_to_num(dset.cable_delay_1))
