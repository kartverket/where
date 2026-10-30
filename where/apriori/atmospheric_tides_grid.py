"""Create interpolators for atmospheric tides coefficients

Description:

Reads atmospheric tides coefficients and creates RectBivariateSpline
interpolators for each dataset.



"""
# External library imports
import numpy as np
from scipy.interpolate import RectBivariateSpline

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where import parsers


@plugins.register
def get_atmospheric_tides(model=None, origin=None):
    """Create interpolators for atmospheric tides coefficients

    Reads the atmospheric tides coefficents from file using the
    AtmosphericTidesParser and creates an RectBiVariateSpline interpolator in
    longtitude and latitude for each coefficient type.

    Longtitude is given from 0 to 360 degrees, while the interpolator should
    work for longtidues from -180 to 180 degrees. The dataset is therefore
    shifted accordingly. The RectBivariateSpline requires longtitude and
    latitude to be in strictly increasing order.

    Returns:
        A dictionary of interpolator functions.
    """
    origin = origin if origin else config.tech.atm_tides.origin.str
    model = model if model else config.tech.atm_tides.model.str
    file_key = "atmospheric_tides_grid_" + model + "_" + origin
    at_data = parsers.parse_key(file_key=file_key).as_dict()

    interpolators = dict()
    lon = np.array(at_data.pop("lon"))
    lat = np.array(at_data.pop("lat"))
    lon = np.unique(lon)
    _, idx = np.unique(lat, return_index=True)

    double_entry = False
    if lon[0] == 0 and lon[-1] == 360:
        # 0 degrees and 360 degrees is the same angle. Remove the 360 data point
        double_entry = True

    # Restore original order
    lat = lat[np.sort(idx)]

    num_value = (len(lat), len(lon))

    # Latitude is given from 90 degrees to -90 degrees
    lat = lat[::-1]
    if double_entry:
        # Strip the last longtitude to avoid double entry for 0 degrees
        lon = lon[:-1]
    # Shift longtitude -180 degrees
    lon = (lon + 180) % 360 - 180
    idx = lon.argsort()
    lon = lon[idx]

    lat = np.radians(lat)
    lon = np.radians(lon)

    for coeff in at_data.keys():
        values = np.array(at_data[coeff]).reshape(num_value)

        # Latitude is given from 90 degrees to -90 degrees
        values = values[::-1, :]
        if double_entry:
            # Strip the last longtitude to avoid double entry for 0 degrees
            values = values[:, :-1]
        # Shift longtitude -180 degrees
        values = values[:, idx]

        interpolators[coeff] = RectBivariateSpline(lon, lat, values.T)

    return interpolators
