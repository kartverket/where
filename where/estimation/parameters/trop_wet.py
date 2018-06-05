"""Calculate the partial derivatives of the wet troposphere

Description:
------------

Calculate the partial derivatives of the troposphere. These will be used to estimate corrections for the wet part of
the tropospheric delay.

References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

.. [2] Kouba, J., Implementation and testing of the gridded Vienna Mapping Function 1 (VMF1),
       Journal of Geodesy(2008) 82:193-205



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
def trop_wet(dset):
    """Calculate the partial derivative of the troposphere for each station

    The partial derivatives are already calculated in the troposphere model and stored in the dataset.

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    stations = np.asarray(dset.unique("station"))

    skip_stations = config.tech.get("skip_stations", section=PARAMETER, default="").list
    skip_idx = [sta == skip_sta for skip_sta in skip_stations for sta in stations]

    if skip_idx:
        stations = stations[np.logical_not(skip_idx)]

    partials = np.zeros((dset.num_obs, len(stations)))
    for multiplier in dset.for_each("station"):
        for site_idx, station in enumerate(stations):
            idx = dset.filter(station=station)
            partials[idx, site_idx] = dset.troposphere_mw[idx] * multiplier

    return partials, stations, "dimensionless"
