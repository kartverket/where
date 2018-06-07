"""Returns correction due to the ionosphere

Description:
------------

The ionospheric correction is already calculated by the correlators and are provided on the NGS file.




"""
# Where imports
from where.lib import plugins


@plugins.register
def ionosphere(dset):
    """Returns the total ionospheric delay for each baseline

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array: Corrections in meters for each observation.
    """
    return dset.iono_delay
