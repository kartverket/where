"""Adds number of scheduled observations to dataset

Description:
------------

Reads the schedule file for the session and adds the number of scheduled observations per station to the dataset

"""

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log
from where import parsers

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def scheduled_obs(dset):
    """Adds scheduled observations to dataset (if available)

    Args:
        dset:     A Dataset containing model data.

    """

    parser = parsers.parse_key("vlbi_schedule_skd", file_vars=config.files.vars)
    
    if parser.data_available:
        data = parser.as_dict()
        log.info(f"{_SECTION}: Adding number of scheduled observations to metadata")
        for sta, num in data.items():
            sta_dict = dset.meta["station"].setdefault(sta, {})
            sta_dict["num_obs_schedule"] = num
    