"""Removes all observations for a given station in a given time interval

Description:
------------

Config should use the following format

    intervals = identifier1 start_epoch1 end_epoch1, identifier2 start_epoch2 end_epoch2, ...

If the identifier is omitted data will be discarded from all stations in the given time interval.
An identifier may be either a station name or a baseline


Example:
--------

    intervals = NYALES20 2013-11-20 17:30:00 2013-11-20 18:24:00
    intervals = 2013-11-20 17:30:00 2013-11-20 18:24:00
    intervals = NYALES20/WETTZELL 2013-11-20 17:30:00 2013-11-20 18:24:00

"""

import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.data.time import Time

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
    scale = config.tech[_SECTION].time_scale.str

    keep_idx = np.ones(dset.num_obs, dtype=bool)
    for interval in intervals:
        interval = interval.split()
        start_time = Time(" ".join(interval[-4:-2]), scale=scale, fmt="iso").datetime
        end_time = Time(" ".join(interval[-2:]), scale=scale, fmt="iso").datetime
        # station name may contain spaces
        station = " ".join(interval[:-4])

        remove_idx = np.logical_and(start_time <= getattr(dset.time, scale).datetime, getattr(dset.time, scale).datetime <= end_time)
        if len(interval) == 5:
            if "/" in station:
                remove_idx = dset.filter(baseline=station, idx=remove_idx)
            else:
                remove_idx = dset.filter(station=station, idx=remove_idx)
        keep_idx = np.logical_and(keep_idx, np.logical_not(remove_idx))

    return keep_idx
