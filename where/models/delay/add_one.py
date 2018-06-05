"""Add one to all corrections

Description:

Silly model just for testing.


$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

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
