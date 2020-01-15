"""Calculate the partial derivatives of the dUT1 Earth Orientation Parameter

Description:
------------

Calculate the partial derivatives of the :math:`UT1 - UTC` Earth orientation parameter.

This is done according to equations (2.34) - (2.36) in Teke :cite:`teke2011`.



"""

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
    column_name = ["dut1"]

    src_dir = dset.src_dir.unit_vector[:, None, :]
    baseline = (dset.site_pos_2.trs.pos - dset.site_pos_1.trs.pos).mat
    dR_dut1 = rotation.dR_dut1(dset.time)
    partials = -(src_dir @ rotation.Q(dset.time) @ dR_dut1 @ rotation.W(dset.time) @ baseline)[:, :, 0]

    return partials, column_name, "meter * (radians per second)"
