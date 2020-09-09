"""Calculate the station-satellite distance rate

Description:
------------

This model computes the rate of the distance between the station and the satellite.

"""
# Standard library
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math import nputil


@plugins.register
def gnss_range_rate(dset: "Dataset"):
    """Calculate rate of the distance between station and satellite in GCRS

    Args:
        rundate:    The model run date.
        tech:       Name of technique.
        dset:       A dataset containing the data.

    Returns:
        table of corrections for each observation
    """
    if "site_vel" not in dset.fields:
        # TODO: This should be replaced by dset.site_posvel
        dset.add_float("site_vel", val=np.zeros([dset.num_obs, 3]), unit="meter/second")
           
    correction =  np.squeeze((dset.sat_posvel.trs.vel - dset.site_vel).val[:, None, :] @ nputil.unit_vector(dset.site_pos.trs.vector_to(dset.sat_posvel.trs.pos))[:, :, None])
    
    return -correction
