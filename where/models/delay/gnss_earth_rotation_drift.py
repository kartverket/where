"""Correct for Earth rotation (Sagnac effect) drift

Description:
------------

This model is used for correcting Doppler shift observations by receiver velocity determinition (in Eq. 6-22 in
 :cite:`zhang2007`). Correction is done either based on precise or broadcast orbits.
"""

# External library imports
import numpy as np

# Midgard imports
from midgard.collections import enums
from midgard.dev import plugins
from midgard.math.constant import constant

# Where imports
from where import apriori


@plugins.register
def gnss_earth_rotation_drift(dset: "Dataset") -> np.ndarray:
    """Determine earth rotation drift based on precise or broadcast satellite clock information

    The correction is caluclated after the description in Section 6.2.9 in :cite:`zhang2007`.

    Args:
        dset (Dataset):   Model data.

    Returns:
        numpy.ndarray:    GNSS earth rotation drift for each observation

    """
    correction = np.zeros(dset.num_obs)

    if "site_vel" not in dset.fields:
        # TODO: This should be replaced by dset.site_posvel
        dset.add_float("site_vel", val=np.zeros([dset.num_obs, 3]), unit="meter/second")

    for sys in dset.unique("system"):
        idx = dset.filter(system=sys)
        omega = constant.get("omega", source=enums.gnss_id_to_reference_system[sys])
        correction[idx] = (
            omega
            / constant.c
            * (
                dset.site_vel[:, 0][idx] * dset.sat_posvel.trs.y[idx]
                - dset.site_vel[:, 1][idx] * dset.sat_posvel.trs.x[idx]
                + dset.site_pos.trs.x[idx] * dset.sat_posvel.trs.vy[idx]
                - dset.site_pos.trs.y[idx] * dset.sat_posvel.trs.vx[idx]
            )
        )

    return correction
