"""Detects outliers based on rms

Description:
------------

"""
import numpy as np

# Where imports
from where.lib import config
from where.lib import plugins

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def rms(dset):
    """Detects outliers based on rms

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    field = config.tech[_SECTION].field.str
    outlier_limit = config.tech[_SECTION].outlier_limit.float
    return np.abs(dset[field]) < outlier_limit * dset.rms(field)
