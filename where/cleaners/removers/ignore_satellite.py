"""Remove all data for given satellites

Description:
------------

Removes all observations of satellites given in the edit file.

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
def ignore_satellite(dset):
    """Edits data based on observing station

    Args:
        dset (Dataset):   A Dataset containing model data.

    Returns:
        numpy.ndarray:    Array containing False for observations to throw away
    """
    satellites = config.tech[_SECTION].satellites.list
    remove_idx = np.zeros(dset.num_obs, dtype=bool)

    if satellites:
        log.info("Discarding observations from satellites: {}", ", ".join(satellites))
        for satellite in satellites:
            remove_idx = np.logical_or(remove_idx, dset.filter(satellite=satellite))

    return np.logical_not(remove_idx)
