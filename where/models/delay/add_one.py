"""Add one to all corrections

Description:

Silly model just for testing.



"""

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import log


@plugins.register
def add_one(dset):
    log.info("Adding one")
    return np.ones(dset.num_obs)
