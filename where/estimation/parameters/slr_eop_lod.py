"""Calculate the partial derivatives of the rate of the Earth Orientation Parameter Length of Day

Description:
------------

Calculate the partial derivatives of the Earth orientation parameter Length of Day.


"""
# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import rotation


@plugins.register
def eop_lod(dset):
    """Calculate the partial derivative of the Earth Orientation Parameter Length of Day

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

    column_name = ["lod"]

    sat_dir = unit_vector[:, None, :]
    dR_dut1 = rotation.dR_dut1(dset.time)

    dt = (dset.time.jd - dset.time.mean.jd)[:, None, None]
    # lod = - ut1_rate * 1 day -> lod_partial = - ut1_rate_partial / 1 day

    pick_relevant_data = np.zeros((dset.num_obs, 3, 3))

    for j in range(0, dset.num_obs):
        if is_core_station[j]:
            pick_relevant_data[j] += np.eye(3)

    partials = (
        sat_dir @ rotation.Q(dset.time) @ dR_dut1 @ rotation.W(dset.time) @ dset.site_pos.trs.val[:, :, None] @ dt
    )[:, :, 0]

    return partials, column_name, "meter * radians / seconds"
