"""Calculate the partial derivatives of the polar motion Earth Orientation Parameters.

Description:
------------

Calculate the partial derivatives of the :math:`x_p` and :math:`y_p` Earth orientation parameters.



"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import rotation


@plugins.register
def eop_pm(dset):
    """Calculate the partial derivative of the polar motion Earth Orientation Parameters

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    # Use only core sites for EOP estimation
    core_sites = config.tech.core_sites.list
    is_core_station = np.in1d(dset.station, core_sites)

    column_names = ["xp", "yp"]
    partials = np.zeros((dset.num_obs, 2))

    unit_vector = dset.sat_pos.gcrs.pos.val - dset.site_pos.gcrs.val
    unit_vector = unit_vector / np.linalg.norm(unit_vector)

    sat_dir = unit_vector[:, None, :]
    pick_relevant_data = np.zeros((dset.num_obs, 3, 3))
    for j in range(0, dset.num_obs):
        if is_core_station[j]:
            pick_relevant_data[j] += np.eye(3)

    # x-pole
    partials[:, 0] = -(
        sat_dir
        @ rotation.Q(dset.time)
        @ rotation.R(dset.time)
        @ rotation.dW_dxp(dset.time)
        @ pick_relevant_data
        @ dset.site_pos.trs.val[:, :, None]
    )[:, 0, 0]

    # y-pole
    partials[:, 1] = -(
        sat_dir
        @ rotation.Q(dset.time)
        @ rotation.R(dset.time)
        @ rotation.dW_dyp(dset.time)
        @ pick_relevant_data
        @ dset.site_pos.trs.val[:, :, None]
    )[:, 0, 0]

    return partials, column_names, "meter per radian"
