"""Clean GNSS observations. 

Description:
------------
Keep only choosen GNSS observations and observation epochs, which at least 4 satellites.

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config

# TODO MURKS: The editor with highest number is processed at the end.
@plugins.register_ordered(
    -99
)  # Before gnss_select_obs.py -> in this case ignore_satellite configuration can be ignored!
def gnss_clean_obs(dset: "Dataset") -> np.ndarray:
    """Clean GNSS observations

    Keep only choosen GNSS observations and observation epochs, which at least 4 satellites.

    NOTE: This only a workaround and works only if observations from one GNSS are used. If several GNSSs are used then
          it can happen that due to the fact that some observation type are only given for a certain GNSS, that all
          observations are deleted.

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away.
    """
    keep_idx = np.ones(dset.num_obs, dtype=bool)

    for epoch in sorted(set(dset.time.gps.mjd)):
        idx = dset.time.gps.mjd == epoch
        num_obs = dset.time.gps.mjd[idx].size
        keep_epoch_idx = np.ones(num_obs, dtype=bool)

        # Keep only satellite observation, which are chosen via configuration
        for obs, sat in enumerate(dset.satellite[idx]):
            if sat[0] not in config.tech.systems.list:
                keep_epoch_idx[obs] = False

        # Reject observation epoch, which includes less than 4 satellite observations
        if np.sum(keep_epoch_idx) < 4:
            keep_epoch_idx = np.zeros(num_obs, dtype=bool)

        # Update final solution
        keep_idx[idx] = keep_epoch_idx

    # +TODO: Workaround -> should be handled by parser or apriori step
    # Loop over GNSSs and observation types
    for sys in dset.meta["obstypes"]:
        for obstype in dset.meta["obstypes"][sys]:

            # Remove observations with (close to) zero value
            keep_idx = np.logical_and(keep_idx, dset.obs[obstype] > 10000)
    # -TODO: Workaround

    return keep_idx
