"""Apply eccentricity vector


Description:
------------

Apply eccentricity vector from the monument to the antenna's reference point (axis intersection)


References:
-----------
@todo - is this described anywhere?



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""

# External library imports
import numpy as np

# WHERE imports
from where import apriori
from where.lib import plugins
from where.lib import log


@plugins.register
def eccentricity_vector(dset):
    """Apply eccentricity vector at all stations.

    Corrections are returned in meters in the Geocentric Celestial Reference System for each
    observation. A Numpy array with 6 columns is returned, the first three columns are \f$ x, y, z \f$ for station 1,
    while the last three columns are \f$ x, y, z \f$ for station 2.

    Args:
        dset:         A Dataset containing model data.

    Returns:
        Numpy array:  GCRS corrections in meters.
    """
    ecc = apriori.get("eccentricity", rundate=dset.rundate)
    data_out = list()
    for _ in dset.for_each("station"):
        data_out.append(eccentricity_vector_station(ecc, dset))

    return np.hstack(data_out)


def eccentricity_vector_station(ecc, dset):
    """Calculate the eccentricity vector for a station.

    Corrections are returned in meters in the Geocentric
    Celestial Reference System for each observation.

    Args:
        dset:        A Dataset containing model data

    Returns:
        Numpy array: GCRS corrections in meters.
    """

    denu = np.full((dset.num_obs, 3), np.nan)
    for site_id in dset.unique("site_id"):
        if ecc[site_id]["type"] == "NEU":
            denu[dset.filter(site_id=site_id)] = ecc[site_id]["vector"][[1, 0, 2]]  # Convert from NEU to ENU

    dxyz = dset.site_pos.convert_enu_to_itrs(denu)
    for site_id in dset.unique("site_id"):
        if ecc[site_id]["type"] == "XYZ":
            dxyz[dset.filter(site_id=site_id)] = ecc[site_id]["vector"]

    missing_index = np.any(np.isnan(dxyz), axis=1)
    for site_id in np.unique(dset.site_id[missing_index]):
        log.error("Missing eccentricity vector for site_id '{}'. Correction set to zero.", site_id)
    dxyz[np.isnan(dxyz)] = 0

    return dset.site_pos.convert_itrs_to_gcrs(dxyz)
