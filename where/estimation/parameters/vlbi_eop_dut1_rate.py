"""Calculate the partial derivatives of the rate of the dUT1 Earth Orientation Parameter

Description:
------------

Calculate the partial derivatives of the rate of the :math:`UT1 - UTC` Earth orientation parameter.


\tau = -\hat{K}QRW\vec{b}


"""

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import rotation


@plugins.register
def eop_dut1_rate(dset):
    """Calculate the partial derivative of the rate of the dUT1 Earth Orientation Parameter

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    column_name = ["ddut1"]

    src_dir = dset.src_dir.unit_vector[:, None, :]
    baseline = (dset.site_pos_2.trs.pos - dset.site_pos_1.trs.pos).mat
    dR_dut1 = rotation.dR_dut1(dset.time)
    dt = (dset.time.jd - dset.time.mean.jd)[:, None, None]
    partials = -(src_dir @ rotation.Q(dset.time) @ dR_dut1 @ rotation.W(dset.time) @ baseline @ dt)[:, :, 0]

    return partials, column_name, "meter * radians * days / seconds"
