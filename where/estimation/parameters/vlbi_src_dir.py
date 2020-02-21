"""Calculate the partial derivatives of the source coordinates

Description:
------------

Calculate the partial derivatives of the source coordinates.

This is done according to equations (2.47) - (2.50) in Teke [2]_.


References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

.. [2] Teke, Kamil, Sub-daily parameter estimation in VLBI data analysis.
       https://geo.tuwien.ac.at/fileadmin/editors/GM/GM87_teke.pdf





"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where import apriori
from where.lib import log

# Name of parameter
PARAMETER = __name__.split(".")[-1]


@plugins.register
def src_dir(dset):
    """Calculate the partial derivative of the source coordinates

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, list of their names, and their unit
    """
    column_names = ["ra", "dec"]
    sources = np.asarray(dset.unique("source"))
    icrf = apriori.get("crf", time=dset.time)

    # Remove sources that should be fixed
    fix_idx = np.zeros(len(sources))

    for group in config.tech[PARAMETER].fix_sources.list:
        fix_idx = np.logical_or(
            [icrf[src].meta[group] if group in icrf[src].meta else src == group for src in sources], fix_idx
        )

    for group in config.tech[PARAMETER].except_sources.list:
        except_idx = np.array([icrf[src].meta[group] if group in icrf[src].meta else src == group for src in sources])
        fix_idx = np.logical_and(np.logical_not(except_idx), fix_idx)

    sources = sources[np.logical_not(fix_idx)]

    # Calculate partials
    partials = np.zeros((dset.num_obs, len(sources) * 2))
    baseline = (dset.site_pos_2.gcrs.pos - dset.site_pos_1.gcrs.pos).mat
    dK_dra = dset.src_dir.dsrc_dra[:, None, :]
    dK_ddec = dset.src_dir.dsrc_ddec[:, None, :]
    all_partials = np.hstack((-dK_dra @ baseline, -dK_ddec @ baseline))[:, :, 0]

    for idx, src in enumerate(sources):
        src_idx = dset.filter(source=src)
        partials[src_idx, idx * 2 : idx * 2 + 2] = all_partials[src_idx]

    column_names = [s + "_" + name for s in sources for name in column_names]

    return partials, column_names, "meter"
