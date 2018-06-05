"""Get apriori data for ocean pole tides

Description:

Read coefficients from file using OceanPoleTidesParser and
construct an interpolator function for each set of
coefficients. Use longtitude and latitude in degrees as data
point coordinates for the interpolator functions.



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""
# External library imports
from scipy.interpolate import RectBivariateSpline
import numpy as np

# Where imports
from where import parsers
from where.lib import plugins


@plugins.register
def get_ocean_pole_tides():
    """Get ocean pole tide coefficients

    Read coefficients from file using OceanPoleTidesParser and construct a
    nearest-neighbour interpolator for each set of coefficients. Use longtitude
    and latitude as data point coordinates for the interpolator functions.

    Longtitude is given from 0 to 360 degrees, while the interpolator should
    work for longtidues from -180 to 180 degrees. The dataset is therefore
    shifted accordingly.

    Returns:
        A dictionary of interpolator functions.
    """
    opt_data = parsers.parse_key(file_key="ocean_pole_tides_cmc").as_dict()

    interpolators = dict()
    lon = np.array(opt_data.pop("lon"))
    lat = np.array(opt_data.pop("lat"))

    lon = np.unique(lon)
    lat = np.unique(lat)

    num_value = (len(lat), len(lon))

    # Shift longtitude -180 degrees
    lon = (lon + 180) % 360 - 180
    idx = lon.argsort()
    lon = lon[idx]

    lat = np.radians(lat)
    lon = np.radians(lon)

    for coeff in opt_data.keys():
        values = np.array(opt_data[coeff]).reshape(num_value)

        # Shift longtitude -180 degrees
        values = values[:, idx]

        interpolators[coeff] = RectBivariateSpline(lon, lat, values.T)

    return interpolators
