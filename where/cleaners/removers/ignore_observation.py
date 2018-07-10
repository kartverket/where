"""Removes a given observation

Description:
------------

One or several observations can be specified in the configuration file, and will be removed from the analysis.

The configuration should use the following format:

    [ignore_observation]
    observations = epoch1 station1, epoch2 station2, ...

In the case of VLBI, specify the baseline instead:

    [ignore_observation]
    observations = epoch1 station1_1/station1_2, epoch2 station2_1/station2_2, ...

If no station is specified, all observations at the given epoch are removed. For VLBI, if only one station is
specified, all observations at the given epoch involving that station are removed.


Example:
--------

    [ignore_observation]
    observations = 2016-12-28 08:49:16 FORTLEZA/MATERA
"""

import numpy as np

# Where imports
from where.lib import config
from where.lib import plugins
from where.lib.time import Time

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def ignore_observation(dset):
    """Removes given observations

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    observations = config.tech[_SECTION].observations.as_list(split_re=", *")
    keep_idx = np.ones(dset.num_obs, dtype=bool)

    for observation in observations:
        date, time, *stations = observation.split()
        epoch = Time(f"{date} {time}", scale="utc", format="iso")
        stations = (" ".join(stations)).split("/")  # station names may contain spaces, split at slash instead

        remove_idx = dset.time.utc.iso == epoch.utc.iso
        for station in stations:
            remove_idx = dset.filter(station=station, idx=remove_idx)
        keep_idx[remove_idx] = False

    return keep_idx
