"""Add one to all corrections

Description:

Silly model just for testing.



"""

# External library imports
import numpy as np

# Where imports
from where.lib import log
from where.lib import plugins


@plugins.register
def add_one(dset):
    log.info("Adding one")
    return np.ones(dset.num_obs)
