"""Framework for detecting outliers

Description:
------------

Each detector should be defined in a separate .py-file. The function inside the .py-file that should be called needs to
be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def rms(dset):
        ...

"""

# External library imports
import numpy as np

# Where imports
from where.lib import config
from where.lib import log
from where.lib import plugins
from where.reports import report


def detect_outliers(config_key, dset):
    """Detect all outliers for a given session

    Args:
        config_key (String):  The configuration key listing which detectors to apply.
        dset (Dataset):       Dataset containing analysis data.
    """
    prefix = config.analysis.get("analysis", default="").str
    log.info(f"Detecting outliers")
    keep_idxs = plugins.call_all(package_name=__name__, config_key=config_key, prefix=prefix, dset=dset)

    all_keep_idx = np.ones(dset.num_obs, dtype=bool)
    for detector, detector_keep_idx in keep_idxs.items():
        log.info(f"Detecting {sum(~detector_keep_idx):5d} outliers based on {detector}")
        report.data("detector_data", dset, detector_name=detector, keep_idx=detector_keep_idx)
        all_keep_idx = np.logical_and(all_keep_idx, detector_keep_idx)

    log.info(f"Removing {sum(~all_keep_idx)} of {dset.num_obs} observations")
    return all_keep_idx
