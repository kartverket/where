"""Remove all data for given satellites

Description:
------------

Removes all observations of satellites given in the edit file.

"""
# Standard library imports
from typing import List, Union

# External library imports
import numpy as np

# Where imports
from where.lib import config
from where.lib import log
from where.lib import plugins

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def ignore_satellite(dset: "Dataset", satellites: Union[List[str], None] = None) -> np.ndarray:
    """Edits data based on observing station

    If 'satellites' argument is given, than it overwrites the 'satellites' option in configuration file.

    Args:
        dset:         A Dataset containing model data.
        satellites:   List with satellites names (e.g. [E01, G09, R10])

    Returns:
        Array containing False for observations to throw away
    """
    satellites = config.tech[_SECTION].satellites.list if satellites is None else satellites
    remove_idx = np.zeros(dset.num_obs, dtype=bool)

    if satellites:
        log.info(f"Discarding observations from satellites: {', '.join(satellites)}")
        for satellite in satellites:
            remove_idx |= dset.filter(satellite=satellite)

    return ~remove_idx
