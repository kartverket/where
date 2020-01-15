"""Framework for detecting outliers

Description:
------------

Each detector should be defined in a separate .py-file. The function inside the .py-file that should be called needs to
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


def apply_outlier_detectors(config_key: str, dset: "Dataset") -> np.ndarray:
    """Apply all outlier detectors for a given session

    Args:
        config_key:  The configuration key listing which detectors to apply.
        dset:        Dataset containing analysis data.
    """
    prefix = dset.vars["pipeline"]
    detectors = config.tech[config_key].list
    log.info(f"Apply outlier detectors")
    keep_idxs = plugins.call_all(package_name=__name__, plugins=detectors, prefix=prefix, dset=dset)

    all_keep_idx = np.ones(dset.num_obs, dtype=bool)
    for detector, detector_keep_idx in keep_idxs.items():
        log.info(f"Detecting {sum(~detector_keep_idx):5d} outliers based on {detector}")
        all_keep_idx = np.logical_and(all_keep_idx, detector_keep_idx)

    log.info(f"Removing {sum(~all_keep_idx)} of {dset.num_obs} observations")
    return all_keep_idx


def apply_outlier_detector(detector: str, dset: "Dataset", **kwargs: Dict[Any, Any]) -> None:
    """Apply defined outlier detector for a given session

    Args:
        detector:  The outlier detector name.
        dset:      Dataset containing analysis data.
        kwargs:    Input arguments to the detector.
    """
    log.info(f"Apply outlier detector {detector!r}")
    keep_idx = plugins.call(package_name=__name__, plugin_name=detector, dset=dset, **kwargs)
    log.info(f"Removing {sum(~keep_idx)} of {dset.num_obs} observations")
    dset.subset(keep_idx)
