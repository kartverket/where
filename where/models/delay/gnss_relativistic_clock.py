"""Correct for relativistic clock effect due to orbit eccentricity

Description:
------------





"""
# Where imports
from where import apriori
from where.lib import log
from where.lib import plugins


@plugins.register
def gnss_relativistic_clock(dset):
    """Determine relativistic clock correction due to orbit eccentricity

    The correction is caluclated for precise and broadcast orbits after Eq. 10.10 and 10.11 in :cite:`iers2010`.

    TODO: 'gnss_relativistic_clock' could be saved in Dataset and used here to avoid calculation of relativistic clock
          correction twice.

    Args:
        dset (Dataset):   Model data.

    Returns:
        numpy.ndarray:    Relativistic clock correction due to orbit eccentricity corrections for each observation

    """
    # TODO: Is there a way that relativistic_clock_correction() could be called directly, without generating an orbit
    #       object? -> Should it be handled over dset.gnss_relativistic_clock ????
    orbit = apriori.get(
        "orbit",
        apriori_orbit=dset.vars["orbit"],
        rundate=dset.rundate,
        time=dset.sat_time,
        satellite=tuple(dset.satellite),
        system=tuple(dset.system),
        station=dset.vars["station"],
    )

    return -orbit.relativistic_clock_correction()
