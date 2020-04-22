"""Framework for removing observations from a Dataset

Description:
------------

Observations that should be discarded can be identified by implementing rejectors. The rejectors can be applied
independently to the Dataset or they can be applied sequentially. In the sequential case the order of the rejectors
matter since observations will be removed immediately after each rejector is finished.

Each rejector should be defined in a separate .py-file. The function inside the .py-file that should be called needs to
be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midgard.dev import plugins

    @plugins.register
    def rms(dset):
        ...

"""
# Standard library imports
from typing import Any, Dict

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log


def apply_observation_rejectors(config_key: str, dset: "Dataset", independent: bool) -> np.ndarray:
    """Apply all configured observation rejectors

    Args:
        config_key:     The configuration key listing which rejectors to apply.
        dset:           Dataset containing analysis data.
        independent:    Flag to indicate whether the rejectors are applied independently or sequentially

    Returns:
        Dataset with rejected observation
    """
    prefix = dset.vars["pipeline"]
    rejectors = config.tech[config_key].list
    word = "independently" if independent else "sequentially"
    num_obs_before = dset.num_obs

    log.info(f"Applying observation rejectors {word}")
    all_keep_idx = np.ones(num_obs_before, dtype=bool)
    for rejector in rejectors:
        rejector_keep_idx = plugins.call(package_name=__name__, plugin_name=rejector, prefix=prefix, dset=dset)
        if independent:
            all_keep_idx = np.logical_and(all_keep_idx, rejector_keep_idx)
        else:
            dset.subset(rejector_keep_idx)
        log.info(f"Found {sum(~rejector_keep_idx):5d} observations based on {rejector}")

    if independent:
        dset.subset(all_keep_idx)
    log.info(f"Removing {num_obs_before - dset.num_obs} of {num_obs_before} observations")
    return dset


def apply_observation_rejector(rejector: str, dset: "Dataset", **kwargs: Dict[Any, Any]) -> None:
    """Apply defined outlier detector for a given session

    Args:
        detector:  The outlier detector name.
        dset:      Dataset containing analysis data.
        kwargs:    Input arguments to the detector.
    """
    log.info(f"Applying observation rejectors {rejector!r}")
    keep_idx = plugins.call(package_name=__name__, plugin_name=rejector, dset=dset, **kwargs)
    log.info(f"Removing {sum(~keep_idx)} of {dset.num_obs} observations")
    dset.subset(keep_idx)
