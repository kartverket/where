"""Apply non tidal atmospheric loading displacements


Description:
------------

"""

# External library imports
import numpy as np

# WHERE imports
from where import apriori
from where.lib import plugins


@plugins.register
def non_tidal_atmospheric_loading(dset):
    """Apply non tidal atmospheric loading displacements at all stations.

    Corrections are returned in meters in the Geocentric Celestial Reference System for each
    observation. A Numpy array with 6 columns is returned, the first three columns are \f$ x, y, z \f$ for station 1,
    while the last three columns are \f$ x, y, z \f$ for station 2.

    Args:
        dset:         A Dataset containing model data.

    Returns:
        Numpy array:  GCRS corrections in meters.
    """
    ntapl = apriori.get("non_tidal_atmospheric_loading", time=dset.time)
    data_out = list()
    for _ in dset.for_each("station"):
        data_out.append(non_tidal_atmospheric_loading_station(ntapl, dset))

    return np.hstack(data_out)


def non_tidal_atmospheric_loading_station(ntapl, dset):
    """Apply non tidal atmospheric loading displacements for a station field.

    Corrections are returned in meters in the Geocentric
    Celestial Reference System for each observation.

    Args:
        dset:        A Dataset containing model data

    Returns:
        Numpy array: GCRS corrections in meters.
    """
    lat, lon, _ = dset.site_pos.llh.T

    dup = ntapl["up"](dset.time, lon, lat)
    deast = ntapl["east"](dset.time, lon, lat)
    dnorth = ntapl["north"](dset.time, lon, lat)

    denu = np.stack((deast, dnorth, dup), axis=1)
    dxyz = dset.site_pos.convert_enu_to_itrs(denu)
    return dset.site_pos.convert_itrs_to_gcrs(dxyz)
