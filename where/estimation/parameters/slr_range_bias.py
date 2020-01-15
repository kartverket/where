"""Calculate the slr station dependent range bias as a constant parameter over the arc

Description:
------------

Calculate the partial derivatives of the measurement with respect to the range bias.

References:
-----------



"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins


# Name of parameter
PARAMETER = __name__.split(".")[-1]


@plugins.register
def slr_range_bias(dset):
    """Calculate the partial derivative of the measurement with respect to range bias.

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple:    Array of partial derivatives, and list of names of derivatives
    """
    # Check if stations to estimate for actually has data
    # this step might be unnecessary, since this is done somwhere else?
    stations = np.asarray(dset.unique("station"))
    stations_to_estimate_for = dset.station[np.asarray(dset.estimate_range, dtype=bool)]
    stations_idx = np.in1d(stations, stations_to_estimate_for)
    stations = stations[stations_idx]

    # Set up the partials of measurement with respect to range bias, which is 1 if station
    # is involved in the measurement, otherwise zero.
    partials = np.zeros((dset.num_obs, len(stations)))

    for site_idx, station in enumerate(stations):
        idx = dset.filter(station=station)
        partials[idx, site_idx] = 1

    column_names = [s for s in stations]
    return partials, column_names, "dimensionless"
