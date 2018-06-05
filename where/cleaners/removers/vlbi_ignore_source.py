"""Remove all data for given sources

Description:
------------

Removes all observations involving sources given in the edit file.

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
def ignore_source(dset):
    """Edits data based on observed source

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    sources = config.tech[_SECTION].sources.list
    remove_idx = np.zeros(dset.num_obs, dtype=bool)

    if sources:
        log.info(f"Discarding observations with sources: {', '.join(sources)}")
        for source in sources:
            remove_idx = np.logical_or(remove_idx, dset.filter(source=source))

    return np.logical_not(remove_idx)
