"""Calculate the station to satellite distance

Description:

For now this is just a SILLY model just for testing.



"""

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins


@plugins.register
def doris_range(dset):
    # Pretend the satellite is "stuck" at a fixed point in space
    sat_pos = np.array([7000000, 0, 0])
    return np.linalg.norm(sat_pos - dset.site_pos.gcrs, axis=1)
