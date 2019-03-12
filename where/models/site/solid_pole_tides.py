"""Calculate the station displacement due to solid pole tides

Description:
------------

Correct for solid pole tide. This model is based on the rotational deformation due to polar motion model as described
by the IERS (see Section 7.1.4 of IERS Conventions[1]_).

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
from where.lib import plugins
from where.lib.unit import Unit


@plugins.register
def solid_pole_tides(dset):
    """Calculate solid pole tide corrections at both stations

    Ocean tide corrections are returned in meters in the Geocentric Celestial Reference System for each observation. A
    Numpy array with 6 columns is returned, the first three columns are \f$ x, y, z \f$ for station 1, while the last
    three columns are \f$ x, y, z \f$ for station 2.

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array with ocean tide corrections in meters.
    """
    data_out = list()
    for _ in dset.for_each("station"):
        data_out.append(solid_pole_tides_station(dset))

    return np.hstack(data_out)


def solid_pole_tides_station(dset):
    """Calculate the solid pole tide corrections for a station

    Solid pole tide corrections are returned in meters in the Geocentric Celestial Reference System for each
    observation.

    Args:
        dset:        A Dataset containing model data.
        version:     Version number of conventional mean pole model

    Returns:
        Numpy array with ocean tide corrections in meters.
    """
    eop = apriori.get("eop", time=dset.time)

    # Equation (7.24) IERS Conventions 2010
    m_1 = eop.x - eop.x_pole
    m_2 = eop.y_pole - eop.y

    lat, lon, _ = dset.site_pos.llh.T
    theta = np.pi / 2 - lat
    coslon, sinlon = np.cos(lon), np.sin(lon)

    # Equation (7.26) IERS Conventions 2010
    dup = -33 * np.sin(2 * theta) * (m_1 * coslon + m_2 * sinlon) * Unit.mm2m
    dsouth = -9 * np.cos(2 * theta) * (m_1 * coslon + m_2 * sinlon) * Unit.mm2m
    deast = 9 * np.cos(theta) * (m_1 * sinlon + m_2 * coslon) * Unit.mm2m

    # Correction in topocentric (east, north, up) coordinates
    denu = np.vstack((deast, -dsouth, dup)).T
    return dset.site_pos.convert_enu_to_gcrs(denu)
