"""Edits data based on elevation

Description:
------------

Identifies observations from the dataset with elevation angle lower than a configured limit.

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def elevation(dset):
    """Edits data based on elevation

    If there are more than one station per observation (i.e. VLBI), all stations need to observe above the threshold
    elevation for the observation to be kept.

    Args:
        dset (Dataset):   A Dataset containing model data.

    Returns:
        numpy.ndarray:    Array containing False for observations to throw away.
    """
    elev_threshold = config.tech[_SECTION].cut_off.float
    keep_idx = np.ones(dset.num_obs, dtype=bool)

    # Avoid rounding errors close to 0 when converting to radians
    if elev_threshold <= 0:
        return keep_idx

    # Check elevation for each station in an observation
    for _ in dset.for_each_suffix("site_pos"):
        keep_idx = np.logical_and(keep_idx, dset.site_pos.elevation >= np.radians(elev_threshold))

    return keep_idx
