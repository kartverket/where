"""Correct for the satellite clock

Description:
------------

Correct satellite clock either based on precise or broadcast orbits.
"""

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori


@plugins.register
def gnss_satellite_clock(dset):
    """Determine satellite clock correction based on precise or broadcast satellite clock information

    Args:
        dset (Dataset):   Model data.

    Returns:
        numpy.ndarray:    GNSS satellite clock corrections for each observation

    """
    if "delay.gnss_satellite_clock" in dset.fields:
        return dset.delay.gnss_satellite_clock
    else:
        orbit = apriori.get(
            "orbit",
            apriori_orbit=dset.vars["orbit"],
            rundate=dset.analysis["rundate"],
            system=tuple(dset.unique("system")),
            station=dset.vars["station"],
        )

        return -orbit.satellite_clock_correction(dset)
