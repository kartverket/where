"""Adds analysis status to dataset

Description:
------------

Adds analysis status to the meta field of the dataset.

"""

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def analysis_status(dset):
    """Adds analysis status to dataset

    Args:
        dset:     A Dataset containing model data.

    """
    status = config.tech[_SECTION].status.str
    dset.meta["analysis_status"] = status
    log.info(f"{_SECTION}: Setting status to '{status}'")
