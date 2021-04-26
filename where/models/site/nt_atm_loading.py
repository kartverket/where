"""Apply non tidal atmospheric loading displacements


Description:
------------

"""

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where.data import position
from where.lib import log

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
    for _ in dset.for_each_suffix("station"):
        data_out.append(non_tidal_atmospheric_loading_station(ntapl, dset))

    return data_out


def non_tidal_atmospheric_loading_station(ntapl, dset):
    """Apply non tidal atmospheric loading displacements for a station field.

    Corrections are returned in meters in the Geocentric
    Celestial Reference System for each observation.

    Args:
        dset:        A Dataset containing model data

    Returns:
        Numpy array: GCRS corrections in meters.
    """
    lat, lon, _ = dset.site_pos.pos.llh.T
    try:
        dup = ntapl["up"](dset.time, lon, lat)
        deast = ntapl["east"](dset.time, lon, lat)
        dnorth = ntapl["north"](dset.time, lon, lat)
    except KeyError:
        log.warn(f"No non-tidal atmospheric loading available for {dset.rundate}")
        dup = np.zeros(len(dset.time))
        deast = np.zeros(len(dset.time))
        dnorth = np.zeros(len(dset.time))

    denu = np.stack((deast, dnorth, dup), axis=1)
    if position.is_position(dset.site_pos):
        pos_correction = position.PositionDelta(denu, system="enu", ref_pos=dset.site_pos, time=dset.time)
    elif position.is_posvel(dset.site_pos):
        # set velocity to zero
        denu = np.concatenate((denu, np.zeros(denu.shape)), axis=1)
        pos_correction = position.PosVelDelta(denu, system="enu", ref_pos=dset.site_pos, time=dset.time)
    else:
        log.fatal(f"dset.site_pos{dset.default_field_suffix} is not a PositionArray or PosVelArray.")

    return pos_correction.gcrs
