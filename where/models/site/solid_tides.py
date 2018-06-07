"""Calculate the station displacement due to solid tides

Description:
------------

Correct for solid tides. IERS Conventions [1]_.

References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
       IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html



"""
# External library imports
import numpy as np

# Where imports
from where import apriori
from where.ext import iers_2010 as iers
from where.lib import plugins


@plugins.register
def solid_tides(dset):
    """Calculate solid tide corrections at both stations

    Solid tide corrections are returned in meters in the Geocentric Celestial Reference System for each
    observation. For each station a Numpy array with 3 columns are created, with :math:`x, y, z`-displacements. Thus,
    for VLBI a Numpy array with 6 columns is returned, the first three columns are :math:`x, y, z` for station 1, while
    the last three columns are :math:`x, y, z` for station 2. For the other techniques, a 3-column Numpy array is
    returned.

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array: Solid tide corrections in meters.

    """
    data_out = list()
    for _ in dset.for_each("station"):
        data_out.append(solid_tides_station(dset))

    return np.hstack(data_out)


def solid_tides_station(dset):
    """Calculate the solid tide corrections for a station

    Solid tide corrections are returned in meters in the Geocentric Celestial Reference System for each observation.

    Args:
        dset:        A Dataset containing model data.

    Returns:
        Numpy array with solid tide corrections in meters.
    """
    eph = apriori.get("ephemerides", time=dset.time)
    dxyz = np.zeros((dset.num_obs, 3))
    obs_dt = dset.time.utc.datetime
    hour_of_day = dset.time.utc.jd_frac * 24

    sun_itrs = eph.pos_itrs("sun")
    moon_itrs = eph.pos_itrs("moon")

    # Calculate correction
    for obs, _ in enumerate(dset.values()):
        dxyz[obs] = iers.dehanttideinel(
            dset.site_pos.itrs_pos[obs],
            obs_dt[obs].year,
            obs_dt[obs].month,
            obs_dt[obs].day,
            hour_of_day[obs],
            sun_itrs[obs],
            moon_itrs[obs],
        )

    # Transform from terrestial to celestial
    return dset.site_pos.convert_itrs_to_gcrs(dxyz)
