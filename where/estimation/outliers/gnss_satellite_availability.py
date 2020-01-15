"""Check GNSS satellite availability

Description:
------------
Check if at least 4 GNSS satellite observation are available in each epoch after outlier rejection.

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins


@plugins.register
def gnss_satellite_availability(dset: "Dataset") -> np.ndarray:
    """Check GNSS satellite availability

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """

    keep_idx = np.ones(dset.num_obs, dtype=bool)
    for epoch in sorted(set(dset.time.gps.mjd)):
        idx = dset.time.gps.mjd == epoch
        num_obs = dset.time.gps.mjd[idx].size
        keep_epoch_idx = np.ones(num_obs, dtype=bool)

        # Reject observation epoch, which includes less than 4 satellite observations
        if np.sum(keep_epoch_idx) < 4:
            keep_epoch_idx = np.zeros(num_obs, dtype=bool)

        # Update final solution
        keep_idx[idx] = keep_epoch_idx

    return keep_idx
