"""Correct for relativistic clock effect due to orbit eccentricity

Description:
------------

Correct relativistic clock effect due to orbit eccentricity either based on precise or broadcast orbits.
"""

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori


@plugins.register
def gnss_relativistic_clock(dset):
    """Determine relativistic clock correction due to orbit eccentricity

    The correction is caluclated for precise and broadcast orbits after Eq. 10.10 and 10.11 in :cite:`iers2010`.

    Args:
        dset (Dataset):   Model data.

    Returns:
        numpy.ndarray:    Relativistic clock correction due to orbit eccentricity corrections for each observation

    """
    if "delay.gnss_relativistic_clock" in dset.fields:
        return dset.delay.gnss_relativistic_clock
    else:
        orbit = apriori.get(
            "orbit",
            apriori_orbit=dset.vars["orbit"],
            rundate=dset.analysis["rundate"],
            system=tuple(dset.unique("system")),
            station=dset.vars["station"],
        )

        return -orbit.relativistic_clock_correction(dset.sat_posvel.trs.pos.val, dset.sat_posvel.trs.vel.val)
