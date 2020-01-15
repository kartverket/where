"""Add dilution of precision to dataset calculated via elevation and azimuth

Description:
------------

Dilution of precision is calculated based on elevation and azimuth between station and satellite for each observation 
epoch. Following dilution of precisions are added:
    GDOP - Geometic DOP
    PDOP - Position DOP
    TDOP - Time DOP
    HDOP - Horizontal DOP
    VDOP - Vertical DOP

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.gnss.compute_dops import compute_dops


# Where imports
from where.lib import log

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def gnss_dop(dset: "Dataset") -> None:
    """Adds dilution of precision (DOP) to dataset

    Args:
        dset:     A Dataset containing model data.
    """
    dops = {
        "gdop": np.zeros((dset.num_obs)),
        "pdop": np.zeros((dset.num_obs)),
        "tdop": np.zeros((dset.num_obs)),
        "hdop": np.zeros((dset.num_obs)),
        "vdop": np.zeros((dset.num_obs)),
    }

    # TODO: Check number of satellite observations !!!
    for time in dset.unique("time"):
        idx = dset.filter(time=time)
        dops["gdop"][idx], dops["pdop"][idx], dops["tdop"][idx], dops["hdop"][idx], dops["vdop"][idx] = compute_dops(
            dset.site_pos.azimuth[idx], dset.site_pos.elevation[idx]
        )

    for dop, val in dops.items():
        if dop in dset.fields:
            dset[dop][:] = val
            log.debug(f"{_SECTION}: Update gdop, pdop, hdop and vdop fields to Dataset.")

        else:
            dset.add_float(dop, val=val)
            log.debug(f"{_SECTION}: Add gdop, pdop, hdop and vdop fields to Dataset.")
