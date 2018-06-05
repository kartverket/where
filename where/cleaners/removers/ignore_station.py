"""Remove all data for given stations

Description:
------------

Removes all observations involving stations given in the edit file.

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
def ignore_station(dset):
    """Edits data based on observing station

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    stations = config.tech[_SECTION].stations.list
    remove_idx = np.zeros(dset.num_obs, dtype=bool)

    if stations:
        log.info("Discarding observations from stations: {}", ", ".join(stations))
        for station in stations:
            remove_idx = np.logical_or(remove_idx, dset.filter(station=station))

    return np.logical_not(remove_idx)
