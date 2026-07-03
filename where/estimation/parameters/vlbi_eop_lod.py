"""Calculate the partial derivatives of the rate of the Earth Orientation Parameter Length of Day

Description:
------------

Calculate the partial derivatives of the Earth orientation parameter Length of Day.





"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import rotation


@plugins.register
def eop_lod(dset):
    """Calculate the partial derivative of the Earth Orientation Parameter Length of Day

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    # Only far field observations is used to estimate this parameter
    idx = ~dset.near_field_obs

    column_name = ["lod"]
    partials = np.zeros((dset.num_obs, 1))

    time = dset.time[idx]
    src_dir = dset.src_dir.unit_vector[:, None, :][idx]
    baseline = (dset.site_pos_2.trs.pos[idx] - dset.site_pos_1.trs.pos[idx]).mat
    dR_dut1 = rotation.dR_dut1(time)
    dt = (time.jd - time.mean.jd)[:, None, None]
    # lod = - ut1_rate * 1 day -> lod_partial = - ut1_rate_partial / 1 day
    partials[idx] = (src_dir @ rotation.Q(time) @ dR_dut1 @ rotation.W(time) @ baseline @ dt)[:, :, 0]

    return partials, column_name, "meter * radians / seconds"
