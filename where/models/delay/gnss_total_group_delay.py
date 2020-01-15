"""Correct for total group delay 

Description:
------------

Total group delay correction is needed for single-frequency positioning using broadcast navigation messages.
"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where.lib import config
from where.lib import log


@plugins.register
def gnss_total_group_delay(dset):
    """Correct for total group delay 

    Args:
        dset (Dataset):   Model data.

    Returns:
        numpy.ndarray:    GNSS total group delay correction for each observation in meter

    """
    if config.tech.apriori_orbit.str != "broadcast":
        log.warn(
            f"Total group delay can only be applied for '{config.tech.apriori_orbit}' orbits. Apriori orbit has "
            f"to be 'broadcast'."
        )
        return np.zeros(dset.num_obs)

    brdc = apriori.get(
        "orbit",
        rundate=dset.analysis["rundate"],
        system=tuple(dset.unique("system")),
        station=dset.vars["station"],
        apriori_orbit="broadcast",
    )

    return brdc.total_group_delay(dset)
