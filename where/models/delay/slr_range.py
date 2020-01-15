"""Calculate the station-satellite distance

Description:
------------

This is the main SLR model, computing the distance between the station and the satellite, taking into account the
gravity field of the Earth.




"""
# Midgard imports
from midgard.dev import plugins

# Where imports
from midgard.math.constant import constant


@plugins.register
def slr_range(dset):
    """Calculate the distance between station and satellite

    Integrate differential equation of motion of the satellite and differential equation of the state transition matrix
    of the satellite.

    Args:
        dset:       A dataset containing the data

    Returns:
        Numpy array: Distance for each observation in meters
    """
    return dset.up_leg * constant.c
