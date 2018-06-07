"""Calculate the station-satellite distance

Description:
------------

This model computes the distance between the station and the satellite.





"""
# External library imports
import numpy as np

# Where imports
from where.lib import plugins


@plugins.register
def gnss_range(dset):
    """Calculate distance between station and satellite in GCRS

    Args:
        rundate:    The model run date.
        tech:       Name of technique.
        dset:       A dataset containing the data.

    Returns:
        table of corrections for each observation
    """

    return np.linalg.norm(dset.sat_posvel.gcrs_pos - dset.site_pos.gcrs, axis=1)
