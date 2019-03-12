"""Detects outliers based on Chi-square test

Description:
------------

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.gnss.solution_validation import sol_validation

# Where imports
from where.lib import config
from where.lib import log
from where.lib import plugins

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def chi2(dset):
    """Detects outliers based on chi-square test

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    field = config.tech[_SECTION].field.str
    alpha = config.tech[_SECTION].alpha.float

    num_params = 4  # TODO
    keep_idx = np.ones(dset.num_obs, dtype=bool)
    for time in dset.unique("time"):
        idx = dset.filter(time=time)
        residual_norm = (dset.residual[idx] - np.mean(dset.residual[idx])) / np.std(dset.residual[idx])
        keep_idx[idx] = sol_validation(residual_norm, alpha, num_params)
    return keep_idx
