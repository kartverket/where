"""Framework for removing observations

Description:
------------

Each remover should be defined in a separate .py-file. The function inside the .py-file that should be called needs to
be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def ignore_station(dset):
        ...

"""

# External library imports
import numpy as np

# Where imports
from where.lib import config
from where.lib import log
from where.lib import plugins
from where.reports import report


def apply_removers(config_key, dset):
    """Apply all removers for a given session

    Args:
        config_key (String):  The configuration key listing which removers to apply.
        dset (Dataset):       Dataset containing analysis data.
    """
    prefix = config.analysis.get("analysis", default="").str
    log.info(f"Applying removers")
    keep_idxs = plugins.call_all(package_name=__name__, config_key=config_key, prefix=prefix, dset=dset)

    all_keep_idx = np.ones(dset.num_obs, dtype=bool)
    for remover, remover_keep_idx in keep_idxs.items():
        log.info(f"Removing {sum(np.logical_not(remover_keep_idx)):5d} observations based on {remover}")
        report.data("remover_data", dset, remover_name=remover, keep_idx=remover_keep_idx)
        all_keep_idx = np.logical_and(all_keep_idx, remover_keep_idx)

    log.info(f"Keeping {sum(all_keep_idx)} of {dset.num_obs} observations")
    dset.subset(all_keep_idx)
