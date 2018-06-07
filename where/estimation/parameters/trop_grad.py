"""Calculate the partial derivatives of horizontal gradients of the troposphere

Description:
------------

Calculate the partial derivatives of the troposphere with respect to the horizontal gradients. These will be used
to estimate corrections of the tropospheric delay.

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
def trop_grad(dset):
    """Calculate the partial derivative of the troposphere with respect to the horizontal gradients for each station

    The components needed to compute the partial derivatives are already calculated in the troposphere model and
    stored in the dataset.

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple:    Array of partial derivatives, and list of names of derivatives
    """
    stations = np.asarray(dset.unique("station"))

    skip_stations = config.tech.get("skip_stations", section=PARAMETER, default="").list
    skip_idx = [sta == skip_sta for skip_sta in skip_stations for sta in stations]

    if skip_idx:
        stations = stations[np.logical_not(skip_idx)]

    partials = np.zeros((dset.num_obs, len(stations) * 2))
    for multiplier in dset.for_each("station"):
        gn = dset.troposphere_mg * np.cos(dset.site_pos.azimuth)
        ge = dset.troposphere_mg * np.sin(dset.site_pos.azimuth)
        for site_idx, station in enumerate(stations):
            idx = dset.filter(station=station)
            partials[idx, site_idx * 2] = gn[idx] * multiplier
            partials[idx, site_idx * 2 + 1] = ge[idx] * multiplier

    column_names = [s + "_" + g for s in stations for g in ["north", "east"]]

    return partials, column_names, "dimensionless"
