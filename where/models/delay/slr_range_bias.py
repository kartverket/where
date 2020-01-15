"""Calculate the station range bias

Description:
------------

This model applies the station range bias
Unit: meters

"""
# Midgard imports
from midgard.dev import plugins


@plugins.register
def slr_range_bias(dset):
    """Calculate the station dependent range bias

    Args:
        dset:       A dataset containing the data

    Returns:
        Numpy array: Bias for each observation in meters
    """
    return dset.range_bias
