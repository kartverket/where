"""Calculate the station to satellite distance

Description:

For now this is just a SILLY model just for testing.


$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

# External library imports
import numpy as np

# Where imports
from where.lib import plugins


@plugins.register
def doris_range(dset):
    # Pretend the satellite is "stuck" at a fixed point in space
    sat_pos = np.array([7000000, 0, 0])
    return np.linalg.norm(sat_pos - dset.site_pos.gcrs, axis=1)
