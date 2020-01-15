"""Apply eccentricity vector


Description:
------------

Apply eccentricity vector from the monument to the antenna's reference point (axis intersection)


References:
-----------
@todo - is this described anywhere?



"""

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where.data import position
from where.lib import log
from where.lib import config

_WARNED_MISSING = set()


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
    ecc = apriori.get("eccentricity", rundate=dset.analysis["rundate"])
    data_out = list()
    for _ in dset.for_each_suffix("station"):
        data_out.append(eccentricity_vector_station(ecc, dset))

    return data_out


def eccentricity_vector_station(ecc, dset):
    """Calculate the eccentricity vector for a station.

    Corrections are returned in meters in the Geocentric
    Celestial Reference System for each observation.

    Args:
        dset:        A Dataset containing model data

    Returns:
        Numpy array: GCRS corrections in meters.
    """
    if position.is_position(dset.site_pos):
        ecc_vector = position.PositionDelta(
            np.zeros((dset.num_obs, 3)), system="enu", ref_pos=dset.site_pos, time=dset.time
        )
    elif position.is_posvel(dset.site_pos):
        ecc_vector = position.PosVelDelta(
            np.zeros((dset.num_obs, 6)), system="enu", ref_pos=dset.site_pos, time=dset.time
        )
    else:
        log.fatal(f"dset.site_pos{dset.default_field_suffix} is not a PositionArray or PosVelArray.")

    fieldnames = config.tech.eccentricity.identifier.list
    fielddata = [dset[field] for field in fieldnames]
    if len(fieldnames) > 1:
        keys = set(tuple(zip(*fielddata)))
    else:
        keys = fielddata[0].tolist()

    for key in keys:
        if len(fieldnames) > 1:
            filters = dict(zip(fieldnames, key))
        else:
            filters = dict(zip(fieldnames, [key]))

        if key not in ecc:
            ecc_vector[dset.filter(**filters), 0:3] = np.zeros(3)
            if key in _WARNED_MISSING:
                continue
            log.warn(f"Missing eccentricity data for {key}. Vector set to zero.")
            _WARNED_MISSING.add(key)
            continue

        if ecc[key]["coord_type"] == "ENU":
            ecc_vector[dset.filter(**filters), 0:3] = ecc[key]["vector"]

    ecc_vector = ecc_vector.trs
    for key in keys:
        if len(fieldnames) > 1:
            filters = dict(zip(fieldnames, key))
        else:
            filters = dict(zip(fieldnames, [key]))

        if key not in ecc:
            ecc_vector[dset.filter(**filters), 0:3] = np.zeros(3)
            continue

        if ecc[key]["coord_type"] == "XYZ":
            ecc_vector[dset.filter(**filters), 0:3] = ecc[key]["vector"]

    return ecc_vector.gcrs
