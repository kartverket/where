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
import math
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.constant import constant


@plugins.register
def slr_relativistic(dset):
    """Calculate relativistic delay for all observations

    Args:
        dset:     Dataset containing the data

    Returns:
        Numpy array:  Correction in meters for each observation
    """
    GM = constant.get("GM", source="egm_2008")
    site_norm = np.linalg.norm(dset.site_pos.gcrs, axis=1)
    sat_norm = np.linalg.norm(dset.sat_pos.gcrs, axis=1)
    sta_sat_norm = np.linalg.norm(dset.sat_pos.gcrs.pos - dset.site_pos.gcrs, axis=1)
    numerator = site_norm + sat_norm + sta_sat_norm
    denominator = site_norm + sat_norm - sta_sat_norm
    correction = np.array([math.log(numerator[i] / denominator[i]) for i in range(0, len(numerator))])

    return 2 * GM / constant.c ** 2 * correction
