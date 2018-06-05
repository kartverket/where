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




$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""
# External library imports
import numpy as np

# Where imports
from where.lib import plugins
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
    icrf = apriori.get("crf", session=dset.dataset_name)

    # Remove sources that should be fixed
    fix_idx = np.zeros(len(sources))

    for group in config.tech[PARAMETER].fix_sources.list:
        fix_idx = np.logical_or(
            [icrf[src].meta[group] if group in icrf[src].meta else src == group for src in sources], fix_idx
        )

    for group in config.tech[PARAMETER].except_sources.list:
        except_idx = (
            np.array([icrf[src].meta[group] if group in icrf[src].meta else src == group for src in sources])
        )
        fix_idx = np.logical_and(np.logical_not(except_idx), fix_idx)

    # Remove sources with few observations
    # TODO redo this type of test after outlier elimination, maybe this should be done in the estimator?
    # Similar test might be useful for stations
    limit = 5
    for idx, src in enumerate(sources):
        src_idx = dset.filter(source=src)
        if np.sum(src_idx) < limit:
            fix_idx[idx] = True
            log.warn("Radio source {} has less than {} observations. Keeping coordinates fixed.".format(src, limit))

    sources = sources[np.logical_not(fix_idx)]

    # Calculate partials
    partials = np.zeros((dset.num_obs, len(sources) * 2))
    zero = np.zeros(dset.num_obs)

    cos_ra = np.cos(dset.src_dir.right_ascension)
    sin_ra = np.sin(dset.src_dir.right_ascension)
    cos_dec = np.cos(dset.src_dir.declination)
    sin_dec = np.sin(dset.src_dir.declination)

    baseline = (dset.site_pos_2.gcrs_pos - dset.site_pos_1.gcrs_pos)[:, :, None]
    dK_dra = np.array([-cos_dec * sin_ra, cos_dec * cos_ra, zero]).T[:, None, :]
    dK_ddec = np.array([-sin_dec * cos_ra, -sin_dec * sin_ra, cos_dec]).T[:, None, :]
    all_partials = np.hstack((dK_dra @ baseline, dK_ddec @ baseline))[:, :, 0]

    for idx, src in enumerate(sources):
        src_idx = dset.filter(source=src)
        partials[src_idx, idx * 2:idx * 2 + 2] = all_partials[src_idx]

    column_names = [s + "_" + name for s in sources for name in column_names]

    return partials, column_names, "meter"
