"""Edits data based on SLR handling file

Description:
------------

Asdf.

"""
# External library imports
import numpy as np

# Where imports
from where import apriori
from where.lib import config
from where.lib import log
from where.lib import plugins


@plugins.register
def data_quality(dset):
    """Edits data based on data quality

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    session = dset.dataset_name
    handling_str = config.session[session].slr_handling.str  # TODO: this is not used ...
    handling = apriori.get("slr_handling_file", time=dset.time)

    remove_idx = np.zeros(dset.num_obs, dtype=bool)
    for station in dset.unique("station"):
        intervals = handling.get(station, {}).get("X", [])
        for interval in intervals:
            start_x = interval[0][0]
            end_x = interval[0][1]
            int_idx = np.logical_and(
                dset.filter(station=station), np.logical_and(dset.date_list >= start_x, dset.date_list <= end_x)
            )
            if np.any(int_idx):
                log.warn("Removed data for station {} in interval {}-{}", station, start_x, end_x)
                remove_idx = np.logical_or(remove_idx, int_idx)
    return np.logical_not(remove_idx)
