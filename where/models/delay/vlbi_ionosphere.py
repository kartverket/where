"""Returns correction due to the ionosphere

Description:
------------

The ionospheric correction may already be calculated and provided in the observation file.




"""
# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import plugins


@plugins.register
def ionosphere(dset):
    """Returns the total ionospheric delay for each baseline

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array: Corrections in meters for each observation.
    """
    # Check if ionosphere is already computed and stored in the dataset
    try:
        return np.nan_to_num(dset.iono_delay)
    except AttributeError:
        pass

    # Check if dTEC is stored on the dataset and try to compute the ionospheric delay
    try:
        print(f"TODO: compute ionosphere")
        return 40.3e16 / np.nan_to_num(dset.ref_freq) ** 2 * np.nan_to_num(dset.dtec)
        # return np.zeros(dset.num_obs)
    except AttributeError:
        # Give up and return zero
        return np.zeros(dset.num_obs)
