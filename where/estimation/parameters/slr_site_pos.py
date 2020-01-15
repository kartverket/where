"""Calculate the partial derivatives of the site positions

Description:
------------

Calculate the partial derivatives of the site positions.


References:
-----------

.. [1] Montenbruck, Oliver and Gill, Eberhard: Satellite Orbits, Springer Verlag, 2000.

.. [2] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html


"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log
from where.lib import rotation
from where.data.time import TimeDelta

# Name of parameter
PARAMETER = __name__.split(".")[-1]


@plugins.register
def site_pos(dset):
    """Calculate the partial derivative of the site position for each station

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple:    Array of partial derivatives, list of their names, and their unit
    """
    # Remove stations that should be fixed
    stations = np.asarray(dset.unique("station"))
    fix_stations = config.tech[PARAMETER].fix_stations.list
    if fix_stations:
        log.info(f"Not estimating stations {fix_stations}")

    # Remove stations with less than 10 observations
    stations_few_obs = []
    for station in stations:
        number_of_observations = np.sum(dset.filter(station=station))
        if number_of_observations < 10:
            stations_few_obs.append(station)

    fix_stations += stations_few_obs
    if stations_few_obs:
        log.info(f"Not estimating stations {stations_few_obs}, less than 10 observations")

    fix_idx = np.in1d(stations, fix_stations)
    if fix_idx.any():
        stations = stations[np.logical_not(fix_idx)]

    reflect_time = dset.time + TimeDelta(dset.time_bias + dset.up_leg, fmt="seconds", scale="utc")
    i2g = rotation.trs2gcrs(reflect_time)
    site_pos_reflect_time = (i2g @ dset.site_pos.val[:, :, None])[:, :, 0]

    unit_vector = dset.sat_pos.pos.val[:] - site_pos_reflect_time[:]
    unit_vector = unit_vector / np.linalg.norm(unit_vector)
    # Calculate partials
    all_partials = -unit_vector[:, None, :] @ i2g
    partials = np.zeros((dset.num_obs, len(stations) * 3))
    for idx, station in enumerate(stations):
        filter_1 = dset.filter(station=station)
        partials[filter_1, idx * 3 : idx * 3 + 3] = all_partials[filter_1][:, 0]

    column_names = [s + "_" + xyz for s in stations for xyz in "xyz"]

    return partials, column_names, "dimensionless"
