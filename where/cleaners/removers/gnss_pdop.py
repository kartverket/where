"""Reject observations based on PDOP

Description:
------------

Identifies observations from the dataset with PDOP greater than a configured limit.

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.cleaners.editors.gnss_dop import gnss_dop
from where.lib import config

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def gnss_pdop(dset: "Dataset") -> np.ndarray:
    """Reject observations based on PDOP

    Args:
        dset:   A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away.
    """
    # Add DOPs to dataset
    gnss_dop(dset)

    # Check PDOP limit
    pdop_limit = config.tech[_SECTION].pdop_limit.float
    keep_idx = np.ones(dset.num_obs, dtype=bool)

    return np.logical_and(keep_idx, dset.pdop < pdop_limit)
