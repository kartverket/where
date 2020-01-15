"""Framework for removing observations

Description:
------------

Each remover should be defined in a separate .py-file. The function inside the .py-file that should be called needs to
be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midgard.dev import plugins

    @plugins.register
    def ignore_station(dset):
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


def apply_removers(config_key: str, dset: "Dataset") -> None:
    """Apply all removers for a given session

    Args:
        config_key:  The configuration key listing which removers to apply.
        dset:        Dataset containing analysis data.
    """
    prefix = dset.vars["pipeline"]
    removers = config.tech[config_key].list
    log.info(f"Applying removers")
    keep_idxs = plugins.call_all(package_name=__name__, plugins=removers, prefix=prefix, dset=dset)

    all_keep_idx = np.ones(dset.num_obs, dtype=bool)
    for remover, remover_keep_idx in keep_idxs.items():
        log.info(f"Removing {sum(np.logical_not(remover_keep_idx)):5d} observations based on {remover}")
        all_keep_idx = np.logical_and(all_keep_idx, remover_keep_idx)

    log.info(f"Keeping {sum(all_keep_idx)} of {dset.num_obs} observations")
    dset.subset(all_keep_idx)

    if dset.num_obs == 0:
        log.fatal("No observations are available.")


def apply_remover(remover: str, dset: "Dataset", **kwargs: Dict[Any, Any]) -> None:
    """Apply defined remover for a given session

    Args:
        remover:   The remover name.
        dset:      Dataset containing analysis data.
        kwargs:    Input arguments to the remover.
    """
    log.info(f"Apply remover {remover!r}")
    keep_idx = plugins.call(package_name=__name__, plugin_name=remover, dset=dset, **kwargs)
    log.info(f"Keeping {sum(keep_idx)} of {dset.num_obs} observations")
    dset.subset(keep_idx)

    if dset.num_obs == 0:
        log.fatal("No observations are available.")
