"""Calculate the partial derivatives of the clock

Description:
------------

Calculate the partial derivatives of the clock. These will be used to estimate corrections for the clock.


References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html




$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""
# External library imports
import numpy as np

# Where imports
from where.lib import plugins
from where.lib import config

# Name of parameter
PARAMETER = __name__.split(".")[-1]


@plugins.register
def clock(dset):
    """Calculate the partial derivative of the clock for each station

    The partial derivatives are simply +/- 1 depending on whether the station is station 1 or station 2 in the
    baseline.

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    stations = dset.unique("station")
    stations.remove(dset.meta["ref_clock"])
    stations = np.asarray(stations)

    skip_stations = config.tech.get("skip_stations", section=PARAMETER, default="").list
    skip_idx = [sta == skip_sta for skip_sta in skip_stations for sta in stations]

    if skip_idx:
        stations = stations[np.logical_not(skip_idx)]

    partials = np.zeros((dset.num_obs, len(stations)))
    for idx, station in enumerate(stations):
        partials[dset.filter(station_1=station), idx] = -1
        partials[dset.filter(station_2=station), idx] = 1

    return partials, stations, "dimensionless"
