"""Adds dilution of precision (DOP) to dataset

Description:
------------

 Add GDOP, PDOP, HDOP and VDOP to dataset.

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.gnss.solution_validation import sol_validation

# Where imports
from where.lib import config
from where.lib import log
from where.lib import plugins

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def gnss_dops(dset: "Dataset") -> None:
    """Adds dilution of precision (DOP) to dataset

    Args:
        dset:     A Dataset containing model data.
    """
    dops = {
        "gdop": np.zeros((dset.num_obs)),
        "pdop": np.zeros((dset.num_obs)),
        "hdop": np.zeros((dset.num_obs)),
        "vdop": np.zeros((dset.num_obs)),
    }

    for time in dset.unique("time"):
        idx = dset.filter(time=time)
        dops["gdop"][idx], dops["pdop"][idx], dops["hdop"][idx], dops["vdop"][idx] = compute_dops(
            dset.site_pos.azimuth[idx], dset.site_pos.elevation[idx]
        )

    for dop, val in dops.items():
        dset.add_float(dop, val=val])

    log.info(f"{_SECTION}: Add gdop, pdop, hdop and vdop fields to Dataset.")
