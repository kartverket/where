"""Calculate the partial derivatives of the vlbi scale

Description:
------------

Calculate the partial derivatives of the vlbi scale.

Implementation is based on cite:`titov2018` equation 10.

"""

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import rotation

# Name of parameter
PARAMETER = __name__.split(".")[-1]


@plugins.register
def scale(dset):
    """Calculate the partial derivative of the vlbi scale

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, list of their names, and their unit
    """
    src_dir = dset.src_dir.unit_vector[:, None, :]
    baseline = (dset.site_pos_2.trs.pos - dset.site_pos_1.trs.pos).mat
    partials = -(src_dir @ rotation.trs2gcrs(dset.time) @ baseline)[:, :, 0]

    return partials, ["scale"], "meter"
