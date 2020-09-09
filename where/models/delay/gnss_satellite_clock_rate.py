"""Correct for satellite clock rate 

Description:
------------

This model is used for correcting Doppler shift observations by receiver velocity determinition (in Eq. 6-14 in
 :cite:`zhang2007`). Correction is done for broadcast orbits. 

TODO: For precise orbits the satellite clock drift and drift rate has to be determined.
"""

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori


@plugins.register
def gnss_satellite_clock_rate(dset):
    """Determine satellite clock rate based on broadcast satellite clock information

    The correction is caluclated after the description in Section 6.2.6.2 in :cite:`zhang2007`.

    Args:
        dset (Dataset):   Model data.

    Returns:
        numpy.ndarray:    GNSS satellite clock rate corrections for each observation

    """
    orbit = apriori.get(
        "orbit",
        apriori_orbit=dset.vars["orbit"],
        rundate=dset.analysis["rundate"],
        system=tuple(dset.unique("system")),
        station=dset.vars["station"],
    )

    return orbit.satellite_clock_rate_correction(dset)
