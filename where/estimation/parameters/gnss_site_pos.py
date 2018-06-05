"""Calculate the partial derivatives of the GNSS site position

Description:
------------

Calculate the partial derivatives of the GNSS site position, e.g. described in :cite:`landau1988` p. 30 and
:cite:`hofmann2008` p. 250.



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""
# External library imports
import numpy as np

# Where imports
from where.lib import plugins


@plugins.register
def gnss_site_pos(dset):
    """Calculate the partial derivative of the GNSS site position

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    station = dset.dataset_name
    partials = np.zeros((dset.num_obs, 3))
    for obs, (sat_pos, site_pos, range_) in enumerate(
        zip(dset.sat_posvel.itrs_pos, dset.site_pos.itrs, dset.gnss_range)
    ):
        partials[obs, 0] = (sat_pos[0] - site_pos[0]) / range_
        partials[obs, 1] = (sat_pos[1] - site_pos[1]) / range_
        partials[obs, 2] = (sat_pos[2] - site_pos[2]) / range_
    column_names = [station + "_x", station + "_y", station + "_z"]

    return partials, column_names, "dimensionless"
