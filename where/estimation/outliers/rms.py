"""Detects outliers based on rms

Description:
------------

"""
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def rms(dset: "Dataset") -> np.ndarray:
    """Detects outliers based on rms

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    field = config.tech[_SECTION].field.str
    outlier_limit = config.tech[_SECTION].outlier_limit.float

    # Epochwise estimation or over whole time period
    if config.tech.estimate_epochwise.bool:
        keep_idx = np.ones(dset.num_obs, dtype=bool)
        for time in dset.unique("time"):
            idx = dset.filter(time=time)
            keep_idx[idx] = np.abs(dset[field][idx]) < outlier_limit * dset.rms(field, idx=idx)
    else:
        keep_idx = np.abs(dset[field]) < outlier_limit * dset.rms(field)

    return keep_idx
