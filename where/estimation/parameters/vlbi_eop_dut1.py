"""Calculate the partial derivatives of the dUT1 Earth Orientation Parameter

Description:
------------

Calculate the partial derivatives of the :math:`UT1 - UTC` Earth orientation parameter.

This is done according to equations (2.34) - (2.36) in Teke :cite:`teke2011`.




$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

# Where imports
from where.ext import sofa_wrapper as sofa
from where.lib import plugins


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
    baseline = (dset.site_pos_2.itrs_pos - dset.site_pos_1.itrs_pos)[:, :, None]
    dR_dut1 = sofa.dR_dut1(dset.time)
    partials = (src_dir @ sofa.Q(dset.time) @ dR_dut1 @ sofa.W(dset.time) @ baseline)[:, :, 0]

    return partials, column_name, "meter * (radians per second)"
