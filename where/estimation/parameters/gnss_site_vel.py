"""Calculate the partial derivatives of the GNSS site velocity

Description:
------------

Calculate the partial derivatives of the GNSS site velocity described in :cite:`petovello2015`.


"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins


@plugins.register
def gnss_site_vel(dset):
    """Calculate the partial derivative of the GNSS site velocity

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    range_ = (dset.sat_posvel.trs.pos - dset.site_pos.trs.pos).length
    partials = -(dset.sat_posvel.pos.trs.val - dset.site_pos.pos.trs.val) / range_[:, None]
    column_names = ["x", "y", "z"]

    return partials, column_names, "dimensionless"
