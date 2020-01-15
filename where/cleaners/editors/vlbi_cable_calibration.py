"""Edit cable calibration data

Description:
------------

Edit cable calibration data based on configuration

"""

# External imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def cable_calibration(dset):
    """Edit cable calibration data

    Args:
        dset:     A Dataset containing model data.

    """
    stations = config.tech[_SECTION].ignore_cable.list

    if stations:
        log.info(f"{_SECTION}: Discarding cable calibration data from {', '.join(stations)}")
        idx_1 = np.zeros(dset.num_obs, dtype=bool)
        idx_2 = np.zeros(dset.num_obs, dtype=bool)

        for station in stations:
            idx_1 = np.logical_or(idx_1, dset.filter(station_1=station))
            idx_2 = np.logical_or(idx_2, dset.filter(station_2=station))

        dset.cable_delay_1[idx_1] = 0.0
        dset.cable_delay_2[idx_2] = 0.0
