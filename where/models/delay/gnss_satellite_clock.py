"""Correct for the satellite clock

Description:
------------

Correct satellite clock either based on precise or broadcast orbits.




$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""
# Where imports
from where import apriori
from where.lib import log
from where.lib import plugins


@plugins.register
def gnss_satellite_clock(dset):
    """Determine satellite clock correction based on precise or broadcast satellite clock information

    TODO: 'gnss_satellite_clock' could be saved in Dataset and used here to avoid calculation of satellite clock
          correction twice.

    Args:
        dset (Dataset):   Model data.

    Returns:
        numpy.ndarray:    GNSS satellite clock corrections for each observation

    """
    orbit = apriori.get(
        "orbit",
        apriori_orbit=dset.vars["orbit"],
        rundate=dset.rundate,
        time=dset.sat_time,
        satellite=tuple(dset.satellite),
        system=tuple(dset.system),
        station=dset.vars["station"],
    )

    return -orbit.satellite_clock_correction()
