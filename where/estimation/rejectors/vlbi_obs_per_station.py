"""Detect observations

Description:
------------

"""
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log
from where.lib import util

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def vlbi_obs_per_station(dset: "Dataset") -> np.ndarray:
    """Detects outliers based on rms

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    keep_idx = np.ones(dset.num_obs, dtype=bool)
    # TODO: simpler test when state fields are constructed better
    if "vlbi_site_pos" not in np.unique(np.char.partition(dset.state.fields, "-")[:, 0]):
        # Keep all observations if station coordinates are not estimated
        return keep_idx

    min_obs = config.tech[_SECTION].min_obs.int
    store_ignore_station = config.tech[_SECTION].store_ignore_station.bool

    stations = dset.unique("station")
    discarded_stations = []
    for station in stations:
        num = dset.num(station=station)
        if num <= min_obs:
            station_idx = dset.filter(station=station)
            keep_idx[station_idx] = False
            discarded_stations.append(station)
    
    if discarded_stations and store_ignore_station:
        log.info(f"Adding {', '.join(discarded_stations)} to ignore_station")
        with config.update_tech_config(dset.analysis["rundate"], dset.vars["pipeline"], session_code=dset.vars["session_code"]) as cfg:
            current = cfg.ignore_station.stations.as_list(", *")
            updated = ", ".join(sorted(current + discarded_stations))
            cfg.update("ignore_station", "stations", updated, source=util.get_program_name())

    return keep_idx
