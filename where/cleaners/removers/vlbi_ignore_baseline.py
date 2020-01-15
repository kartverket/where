"""Remove all observations for given baselines

Description:
------------

Removes all observations for baselines given in the edit file.

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def ignore_baseline(dset):
    """Edits data based on baselines

    Args:
        dset (Dataset):   A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    baselines = config.tech[_SECTION].baselines.as_list(split_re=", *")
    remove_idx = np.zeros(dset.num_obs, dtype=bool)

    if baselines:
        log.info(f"Discarding observations with baselines: {', '.join(baselines)}")

        # Add baselines with stations in reverse order
        baselines.extend(["/".join(reversed(b.split("/"))) for b in baselines])
        for baseline in baselines:
            remove_idx = np.logical_or(remove_idx, dset.filter(baseline=baseline))

    return np.logical_not(remove_idx)
