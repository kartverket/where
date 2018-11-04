"""Returns delay correction due to the cables

Description:
------------

The cable calibration is already calculated by the correlators and is provided on the NGS file.

"""
# External library imports
import numpy as np

# Where imports
from where.lib import config
from where.lib import log
from where.lib import plugins

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
    stations = config.tech[_SECTION].ignore_cable.list

    if stations:
        log.info("Discarding cable calibration data from {}", ", ".join(stations))
        idx_1 = np.zeros(dset.num_obs, dtype=bool)
        idx_2 = np.zeros(dset.num_obs, dtype=bool)

        for station in stations:
            idx_1 = np.logical_or(idx_1, dset.filter(station_1=station))
            idx_2 = np.logical_or(idx_2, dset.filter(station_2=station))

        dset.cable_delay_1[idx_1] = 0.0
        dset.cable_delay_2[idx_2] = 0.0

    return -(np.nan_to_num(dset.cable_delay_2) - np.nan_to_num(dset.cable_delay_1))
