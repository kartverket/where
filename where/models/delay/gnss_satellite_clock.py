"""Correct for the satellite clock

Description:
------------

Correct satellite clock either based on precise or broadcast orbits.
"""
# Where imports
from where import apriori
from where.lib import log
from where.lib import plugins


@plugins.register
def gnss_satellite_clock(dset):
    """Determine satellite clock correction based on precise or broadcast satellite clock information

    Args:
        dset (Dataset):   Model data.

    Returns:
        numpy.ndarray:    GNSS satellite clock corrections for each observation

    """
    if "gnss_satellite_clock" in dset.fields:
        return -dset.gnss_satellite_clock
    else:
        orbit = apriori.get(
            "orbit",
            apriori_orbit=dset.vars["orbit"],
            rundate=dset.rundate,
            time=dset.sat_time,
            satellite=tuple(dset.satellite),
            system=tuple(dset.system),
            station=dset.vars["station"],
        )

        return -orbit.satellite_clock_correction(dset)
