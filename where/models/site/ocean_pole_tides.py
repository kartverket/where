"""Calculate the station displacement due to ocean pole tides

Description:
------------

Correct for ocean pole tide. This model is based on the rotational deformation due to polar motion model as described
by the IERS (see Section 7.1.5 of IERS Conventions [1]_).

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
from where.lib import constant
from where.lib import plugins
from where.lib.unit import unit


@plugins.register
def ocean_pole_tides(dset):
    """Calculate ocean pole tide corrections at all stations

    Ocean pole tide corrections are returned in meters in the Geocentric Celestial Reference System for each
    observation. For each station a Numpy array with 3 columns are created, with :math:`x, y, z`-displacements. Thus,
    for VLBI, a Numpy array with 6 columns is returned, the first three columns are :math:`x, y, z`e for station 1,
    while the last three columns are :math:`x, y, z` for station 2. For other techniques, a 3-column Numpy array is
    returned.

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array: Ocean pole tide corrections in meters.

    """
    opt = apriori.get("ocean_pole_tides")

    data_out = list()
    for _ in dset.for_each("station"):
        data_out.append(ocean_pole_tides_station(dset, opt))

    return np.hstack(data_out)


def ocean_pole_tides_station(dset, opt):
    """Calculate the ocean pole tide corrections for a station

    Ocean pole tide corrections are returned in meters in the Geocentric Celestial Reference System for each
    observation.

    Args:
        dset:        A Dataset containing model data.
        version:     Version number of conventional mean pole model
        opt:         Apriori ocean pole tide coefficients

    Returns:
        Numpy array with ocean tide corrections in meters.
    """
    eop = apriori.get("eop", time=dset.time)

    # Constants in IERS Conventions equation 7.29
    H_p = np.sqrt(8 * np.pi / 15) * constant.omega ** 2 * constant.a ** 4 / constant.GM
    K = 4 * np.pi * constant.G * constant.a * constant.rho_w * H_p / (3 * constant.g_E)
    # @todo make these global? related to love numbers somehow
    gamma_2_R = 0.6870
    gamma_2_I = 0.0036

    # Equation (7.24) IERS Conventions 2010
    m_1 = (eop.x - eop.x_mean) * unit.arcsec2rad
    m_2 = (eop.y_mean - eop.y) * unit.arcsec2rad

    lat, lon, _ = dset.site_pos.llh.T

    u_enu_R = np.array(
        [opt["u_e_R"](lon, lat, grid=False), opt["u_n_R"](lon, lat, grid=False), opt["u_r_R"](lon, lat, grid=False)]
    )
    u_enu_I = np.array(
        [opt["u_e_I"](lon, lat, grid=False), opt["u_n_I"](lon, lat, grid=False), opt["u_r_I"](lon, lat, grid=False)]
    )

    # Equation (7.29) IERS Conventions 2010
    denu = K * ((m_1 * gamma_2_R + m_2 * gamma_2_I) * u_enu_R + (m_1 * gamma_2_R - m_2 * gamma_2_I) * u_enu_I)

    return dset.site_pos.convert_enu_to_gcrs(denu.T)
