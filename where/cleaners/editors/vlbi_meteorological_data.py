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
def meteorological_data(dset):
    """Edit cable calibration data

    Args:
        dset:     A Dataset containing model data.

    """
    pressure_stations = config.tech[_SECTION].ignore_pressure.list
    temperature_stations = config.tech[_SECTION].ignore_temperature.list

    if pressure_stations:
        log.info(f"{_SECTION}: Discarding pressure data from {', '.join(pressure_stations)}")
        idx_1 = np.zeros(dset.num_obs, dtype=bool)
        idx_2 = np.zeros(dset.num_obs, dtype=bool)

        for station in pressure_stations:
            idx_1 = np.logical_or(idx_1, dset.filter(station_1=station))
            idx_2 = np.logical_or(idx_2, dset.filter(station_2=station))

        dset.pressure_1[idx_1] = np.nan
        dset.pressure_2[idx_2] = np.nan

    if temperature_stations:
        log.info(f"{_SECTION}: Discarding temperature data from {', '.join(temperature_stations)}")
        idx_1 = np.zeros(dset.num_obs, dtype=bool)
        idx_2 = np.zeros(dset.num_obs, dtype=bool)

        for station in temperature_stations:
            idx_1 = np.logical_or(idx_1, dset.filter(station_1=station))
            idx_2 = np.logical_or(idx_2, dset.filter(station_2=station))

        dset.temperature_1[idx_1] = np.nan
        dset.temperature_2[idx_2] = np.nan
