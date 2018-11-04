"""Edits data based on data quality

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
def data_quality(dset):
    """Edits data based on data quality

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    dq_threshold = config.tech[_SECTION].threshold.int
    return np.nan_to_num(dset.data_quality) <= dq_threshold
