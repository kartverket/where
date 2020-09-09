"""Correct for relativistic clock rate effect due to orbit eccentricity

Description:
------------

Correct relativistic clock rate effect due to orbit eccentricity based on broadcast orbits.

TODO: for precise orbit the semi-major axis of the orbit has to be determined. Conversion from xyz vxvyvz to kepler.
"""

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.constant import constant

# Where imports
from where import apriori

@plugins.register
def gnss_relativistic_clock_rate(dset):
    """Determine relativistic clock rate correction due to orbit eccentricity

    The correction is caluclated for precise and broadcast orbits after Eq. 6-19 in :cite:`zhang2007`.

    Args:
        dset (Dataset):   Model data.

    Returns:
        numpy.ndarray:    Relativistic clock rate correction due to orbit eccentricity corrections for each observation

    """
    
    gm = constant.get("GM", source="iers_2010")
    
    orbit = apriori.get(
        "orbit",
        apriori_orbit=dset.vars["orbit"],
        rundate=dset.analysis["rundate"],
        system=tuple(dset.unique("system")),
        station=dset.vars["station"],
    )
    
    a = np.power(orbit.get_ephemeris_field("sqrt_a", dset), 2)
    
    correction = 2 * gm/constant.c * (1/a - 1/dset.sat_posvel.trs.pos.length)

    return correction