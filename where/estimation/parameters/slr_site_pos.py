"""Calculate the partial derivatives of the site positions

Description:
------------

Calculate the partial derivatives of the site positions.


References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html


"""
# External library imports
import numpy as np

# Where imports
from where.lib import plugins
from where.lib import config

# Name of parameter
PARAMETER = __name__.split(".")[-1]


@plugins.register
def site_pos(dset):
    """Calculate the partial derivative of the site position for each station

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, list of their names, and their unit
    """
    # Remove stations that should be fixed
    stations = np.asarray(dset.unique("station"))
    fix_stations = config.tech[PARAMETER].fix_stations.list
    fix_idx = np.in1d(stations, fix_stations)
    if fix_idx.any():
        stations = stations[np.logical_not(fix_idx)]

    unit_vector = dset.sat_pos.gcrs[:] - dset.site_pos.gcrs[:]
    # Calculate partials
    all_partials = -unit_vector[:, None, :] @ dset.time.itrs2gcrs
    partials = np.zeros((dset.num_obs, len(stations) * 3))
    for idx, station in enumerate(stations):
        filter_1 = dset.filter(station=station)
        partials[filter_1, idx * 3 : idx * 3 + 3] = all_partials[filter_1][:, 0] * -1

    column_names = [s + "_" + xyz for s in stations for xyz in "xyz"]

    return partials, column_names, "dimensionless"
