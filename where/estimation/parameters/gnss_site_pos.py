"""Calculate the partial derivatives of the GNSS site position

Description:
------------

Calculate the partial derivatives of the GNSS site position, e.g. described in :cite:`landau1988` p. 30 and
:cite:`hofmann2008` p. 250.




"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins


@plugins.register
def gnss_site_pos(dset):
    """Calculate the partial derivative of the GNSS site position

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    partials = (dset.sat_posvel.pos.trs.val - dset.site_pos.pos.trs.val) / dset.delay.gnss_range[:, None]
    column_names = ["x", "y", "z"]

    return partials, column_names, "dimensionless"
