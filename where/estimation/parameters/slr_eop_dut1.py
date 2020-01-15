"""Calculate the partial derivatives of the dUT1 Earth Orientation Parameter

Description:
------------

Calculate the partial derivatives of the :math:`UT1 - UTC` Earth orientation parameter.

This is done according to equations (2.34) - (2.36) in Teke :cite:`teke2011`.



"""
# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import rotation


@plugins.register
def eop_dut1(dset):
    """Calculate the partial derivative of the dUT1 Earth Orientation Parameter

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """

    unit_vector = dset.sat_pos.gcrs.pos.val - dset.site_pos.gcrs.val
    unit_vector = unit_vector / np.linalg.norm(unit_vector)

    # Use only core sites for EOP estimation
    core_sites = config.tech.core_sites.list
    is_core_station = np.in1d(dset.station, core_sites)

    column_name = ["dut1"]

    partials = np.zeros(dset.num_obs)

    sat_dir = unit_vector[:, None, :]
    dR_dut1 = rotation.dR_dut1(dset.time)

    pick_relevant_data = np.zeros((dset.num_obs, 3, 3))

    for j in range(0, dset.num_obs):
        if is_core_station[j]:
            pick_relevant_data[j] += np.eye(3)

    partials = (
        sat_dir
        @ rotation.Q(dset.time)
        @ dR_dut1
        @ rotation.W(dset.time)
        @ pick_relevant_data
        @ dset.site_pos.trs.val[:, :, None]
    )[:, :, 0]

    return partials, column_name, "meter * (radians per second)"
