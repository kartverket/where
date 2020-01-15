"""Edits data based on SLR handling file

Description:
------------



"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where.lib import log


@plugins.register
def data_handling(dset):
    """Edits data based on SLR handling file

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    handling = apriori.get("slr_handling_file", time=dset.time)

    remove_idx = np.zeros(dset.num_obs, dtype=bool)
    for station in dset.unique("station"):
        # TODO: To be implemented
        if "V" in handling.get(station, {}):
            log.dev(f"TODO: Station {station}, marked with a V, not sure what that means")

        # X is data to be deleted
        # N is a non reliable station, not to be used for operational analysis
        # Q is a station in quarantene
        for key in ["X", "N", "Q"]:
            intervals = handling.get(station, {}).get(key, [])
            for interval in intervals:
                start_x, end_x = interval[0]
                int_idx = (
                    dset.filter(station=station) & (dset.time.datetime >= start_x) & (dset.time.datetime <= end_x)
                )
                if np.any(int_idx):
                    log.debug(f"Removed data for station {station} in interval {start_x}-{end_x}, marked with {key}")
                    remove_idx |= int_idx
    return ~remove_idx
