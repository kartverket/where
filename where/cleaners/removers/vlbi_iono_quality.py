"""Edits data based on quality flag of ionosphere correction

Description:
------------

"""
# Where imports
from where.lib import config
from where.lib import plugins

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def iono_quality(dset):
    """Edits data based on iono quality

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    iono_threshold = config.tech[_SECTION].threshold.int
    return dset.iono_quality >= iono_threshold
