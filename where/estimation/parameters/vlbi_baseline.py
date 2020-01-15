"""Calculate the partial derivatives of the vlbi baselines

Description:
------------

Calculate the partial derivatives of the vlbi baselines.


References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html
"""

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Name of parameter
PARAMETER = __name__.split(".")[-1]


@plugins.register
def baseline(dset):
    """Calculate the partial derivative of the vlbi baselines

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, list of their names, and their unit
    """
    baselines = np.asarray(dset.unique("baseline"))

    # Calculate partials
    bs = (dset.site_pos_2 - dset.site_pos_1).pos
    all_partials = -dset.src_dir.trs.val[:, None, :] @ (bs.val / bs.length[:, None])[:, :, None]

    partials = np.zeros((dset.num_obs, len(baselines)))
    for bs_idx, b in enumerate(baselines):
        dset_idx = dset.filter(baseline=b)
        partials[dset_idx, bs_idx] = all_partials[dset_idx, 0, 0]

    return partials, baselines, "dimensionless"
