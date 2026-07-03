"""Calculate the partial derivatives of the dUT1 Earth Orientation Parameter

Description:
------------

Calculate the partial derivatives of the :math:`UT1 - UTC` Earth orientation parameter.

This is done according to equations (2.34) - (2.36) in Teke :cite:`teke2011`.



"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import rotation


@plugins.register
def eop_dut1(dset):
    """Calculate the partial derivative of the dUT1 Earth Orientation Parameter

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    # Only far field observations is used to estimate this parameter
    idx = ~dset.near_field_obs

    column_name = ["dut1"]
    partials = np.zeros((dset.num_obs, 1))

    time = dset.time[idx]
    src_dir = dset.src_dir.unit_vector[:, None, :][idx]
    baseline = (dset.site_pos_2.trs.pos[idx] - dset.site_pos_1.trs.pos[idx]).mat
    dR_dut1 = rotation.dR_dut1(time)
    partials[idx] = -(src_dir @ rotation.Q(time) @ dR_dut1 @ rotation.W(time) @ baseline)[:, :, 0]

    return partials, column_name, "meter * (radians per second)"
