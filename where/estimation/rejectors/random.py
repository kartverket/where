"""Discards a random percentage of observations. 

For testing purposes.

Description:
------------

"""
import random as random_

import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def random(dset: "Dataset") -> np.ndarray:
    """Selects outliers randomly

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    percent = config.tech[_SECTION].percent.float

    # Epochwise estimation or over whole time period
    keep_idx = np.ones(dset.num_obs, dtype=bool)
    discard_num = int(dset.num_obs*percent/100)
    discard_idx = random_.sample(range(dset.num_obs), discard_num)
    keep_idx[discard_idx] = False 
    return keep_idx
