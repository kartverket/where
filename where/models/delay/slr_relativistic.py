"""Calculate relativistic delay

Description:
------------

This model determines the relativistic correction in the laser ranging, see equation (11.17) in [1].
Only the gravity of Earth is considered.

References:
-----------

[1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
    IERS Technical Note No. 36, BKG (2010)



"""
# External library imports
import numpy as np
import math

# Where imports
from where.lib import constant
from where.lib import plugins
from where.lib.time import TimeDelta
from where import apriori


@plugins.register
def slr_relativistic(dset):
    """Calculate relativistic delay for all observations

    Args:
        dset:     Dataset containing the data

    Returns:
        Numpy array:  Correction in meters for each observation
    """
    eph = apriori.get("ephemerides")
    GM = constant.get("GM", source="egm_2008")
    site_norm = np.linalg.norm(dset.site_pos.gcrs, axis=1)
    sat_norm = np.linalg.norm(dset.sat_pos.gcrs, axis=1)
    geocenter_movement = eph.pos_bcrs("earth", time=dset.time + TimeDelta(dset.up_leg, format="sec")) - eph.pos_bcrs(
        "earth", time=dset.time
    )
    sta_sat_norm = np.linalg.norm(dset.sat_pos.gcrs - dset.site_pos.gcrs + geocenter_movement, axis=1)
    numerator = site_norm + sat_norm + sta_sat_norm
    denominator = site_norm + sat_norm - sta_sat_norm
    correction = np.array([math.log(numerator[i] / denominator[i]) for i in range(0, len(numerator))])

    return 2 * GM / constant.c ** 2 * correction
