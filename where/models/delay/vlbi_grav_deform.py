"""Calculate the delay caused by the gravitational deformation

Description:
------------


References:
-----------

"""
# Standard library import

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori


@plugins.register
def vlbi_gravitational_deformation(dset):
    """Calculate gravitational deformation at both stations

    Args:
        dset (Dataset): Model data.

    Returns:
        Numpy array: Corrections in meters for each observation
    """
    data_out = np.zeros(dset.num_obs)
    for multiplier in dset.for_each_suffix("station"):
        data_out += multiplier * gravitational_deformation_station(dset)

    return data_out


def gravitational_deformation_station(dset):
    """Calculate gravitational deformation at one station

    Args:
        dset:        A Dataset containing model data.

    Returns:
        Numpy array: delay due to gravitational deformation in meters.
    """
    deform = apriori.get("vlbi_gravitational_deformation", rundate=dset.analysis["rundate"])
    delays = np.zeros(dset.num_obs)

    for station in dset.unique("station"):
        if station not in deform:
            continue

        sta_idx = dset.filter(station=station)
        delays[sta_idx] = deform[station](dset.site_pos.elevation[sta_idx])
    return delays
