"""Returns corrections in the geometric delay due to the propagation through the atmosphere.

Description:
------------

Calculate the geometric propagation delay using the Consensus model as described in the IERS Conventions
:cite:`iers2010`, section 11.1.




"""
# Midgard imports
from midgard.dev import plugins

# Where imports
from midgard.math.constant import constant
from where.lib import log


@plugins.register_ordered(1000)
def geometric_delay(dset):
    """Returns the part of the geometric delay due to propagation through the atmosphere for each baseline

    This model depends on `troposphere_radio` already having run. Thus, the sort value is set to 1000 to make sure it
    runs last.

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array: Corrections in meters for each observation.

    """
    # Geometric delay due to the atmosphere in equation (11.11)
    if "troposphere_dT_1" in dset.fields:
        datm1 = dset.troposphere_dT_1
    else:
        log.warn("Missing troposphere data. Correction set to zero")
        datm1 = 0

    baseline_gcrs_vel = (dset.site_pos_2.gcrs - dset.site_pos_1.gcrs).vel
    return datm1 * (baseline_gcrs_vel.val[:, None, :] @ dset.src_dir.unit_vector[:, :, None] / constant.c)[:, 0, 0]
