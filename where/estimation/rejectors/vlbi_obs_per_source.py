"""Detect observations

Description:
------------

"""
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def vlbi_obs_per_source(dset: "Dataset") -> np.ndarray:
    """Detects outliers based on rms

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    keep_idx = np.ones(dset.num_obs, dtype=bool)
    # TODO: simpler test when state fields are constructed better
    if "vlbi_src_dir" not in np.unique(np.char.partition(dset.state.fields, "-")[:, 0]):
        # Keep all observations if source coordinates are not estimated
        return keep_idx

    sources = dset.unique("source")
    for source in sources:
        num = dset.num(source=source)
        # At least four good observations is need to determine source positions (two unknowns)
        if num <= 4:
            source_idx = dset.filter(source=source)
            keep_idx[source_idx] = False
    return keep_idx
