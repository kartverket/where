"""Removes observations that are NotANumber

Description:
------------

"""
import numpy as np

# Midgard imports
from midgard.dev import plugins


@plugins.register
def vlbi_ignore_nan(dset):
    """Edits data based on nan

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    return np.logical_not(np.isnan(dset.observed_delay))
