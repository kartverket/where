"""Remove all data for given GNSS

Description:
------------

Removes all observations of GNSS given in the observation data file.

"""
# Standard library imports
from typing import List, Union

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def gnss_ignore_system(dset: "Dataset", systems: Union[List[str], None] = None) -> np.ndarray:
    """Edits data based on observing station

    Args:
        dset:       A Dataset containing model data.
        systems:    List with GNSS identifier (e.g. [G, E])

    Returns:
        Array containing False for observations to throw away
    """
    systems = config.tech[_SECTION].systems.list if systems is None else systems
    remove_idx = np.zeros(dset.num_obs, dtype=bool)

    if systems:
        log.info(f"Discarding observations from GNSS: {', '.join(systems)}")
        for system in systems:
            remove_idx |= dset.filter(system=system)

    # Remove unneccessary meta entries
    if "R" in systems:
        if "glonass_slot" in dset.meta.keys():
            del dset.meta["glonass_slot"]

        if "glonass_bias" in dset.meta.keys():
            del dset.meta["glonass_bias"]

    for sys in systems:
        if sys in dset.meta["obstypes"]:
            del dset.meta["obstypes"][sys]
        if sys in dset.meta["phase_shift"]:
            del dset.meta["phase_shift"][sys]

    return ~remove_idx
