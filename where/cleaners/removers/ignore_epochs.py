"""Removes all observations for a given station in a given time interval

Description:
------------

Config should use the following format

    ignore_epochs = station1 start_epoch1 end_epoch1, station1 start_epoch2 end_epoch2, ...

If the station name is omitted data will be discarded from all stations in the given time interval


Example:
--------

    ignore_epochs = NYALES20 2013-11-20 17:30:00 2013-11-20 18:24:00

"""

import numpy as np

# Where imports
from where.lib import config
from where.lib import plugins
from where.lib.time import Time

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def ignore_epochs(dset):
    """Edits data based on data quality

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    intervals = config.tech[_SECTION].intervals.as_list(split_re=", *")

    keep_idx = np.ones(dset.num_obs, dtype=bool)
    for interval in intervals:
        interval = interval.split()
        start_time = Time(" ".join(interval[-4:-2]), scale="utc", format="iso")
        end_time = Time(" ".join(interval[-2:]), scale="utc", format="iso")
        # station name may contain spaces
        station = " ".join(interval[:-4])

        remove_idx = np.logical_and(start_time < dset.time, dset.time < end_time)
        if len(interval) == 5:
            remove_idx = dset.filter(station=station, idx=remove_idx)
        keep_idx = np.logical_and(keep_idx, np.logical_not(remove_idx))

    return keep_idx
