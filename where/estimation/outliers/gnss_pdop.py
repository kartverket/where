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
from where.cleaners.removers.gnss_pdop import gnss_pdop as gnss_pdop_
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
    return gnss_pdop_(dset)
