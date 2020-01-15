"""Calculate the partial derivatives of the rate of the Earth Orientation Parameter Length of Day

Description:
------------

Calculate the partial derivatives of the Earth orientation parameter Length of Day.





"""

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
    column_name = ["lod"]

    src_dir = dset.src_dir.unit_vector[:, None, :]
    baseline = (dset.site_pos_2.trs.pos - dset.site_pos_1.trs.pos).mat
    dR_dut1 = rotation.dR_dut1(dset.time)
    dt = (dset.time.jd - dset.time.mean.jd)[:, None, None]
    # lod = - ut1_rate * 1 day -> lod_partial = - ut1_rate_partial / 1 day
    partials = (src_dir @ rotation.Q(dset.time) @ dR_dut1 @ rotation.W(dset.time) @ baseline @ dt)[:, :, 0]

    return partials, column_name, "meter * radians / seconds"
