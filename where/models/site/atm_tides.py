"""Calculate the station displacement due to atmospheric tides

Description:
------------

Correct for atmospheric tides. This model is based on the S1 and S2 atmospheric pressure loading described by the IERS
(see Section 7.1.3 of IERS Conventions [1]_).

References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
       IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

.. [2] http://geophy.uni.lu/ggfc-atmosphere/tide-loading-calculator.html




"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where import apriori
from where.lib import config
from where.lib import log
from where.data import position


@plugins.register
def atmospheric_tides(dset):
    """Calculate atmospheric tide corrections at both stations

    Atmospheric tide corrections are returned in meters in the Geocentric Celestial Reference System for each
    observation. A Numpy array with 6 columns is returned, the first three columns are \f$ x, y, z \f$ for station 1,
    while the last three columns are \f$ x, y, z \f$ for station 2.

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array: Atmospheric tide corrections in meters.
    """
    data_out = list()
    for _ in dset.for_each_suffix("station"):
        data_out.append(atmospheric_tides_station(dset))

    return data_out


def atmospheric_tides_station(dset):
    """Calculate the atmospheric tides corrections for a station

    Atmospheric tides corrections are returned in meters in the Geocentric Celestial Reference System for each
    observation.

    Args:
        dset:        A Dataset containing model data

    Returns:
        Numpy array with atmospheric tide corrections in meters.

    """
    coeff = apriori.get("atmospheric_tides")
    use_cmc = config.tech.atmospheric_tides_cmc.bool

    # S1 has a period of 1 cycle/day, S2 has a period of 2 cycle/day
    omega_1 = 2 * np.pi
    omega_2 = 4 * np.pi

    # Time argument is fraction of UT1 day, see [2].
    t = dset.time.ut1.jd_frac
    lat, lon, _ = dset.site_pos.pos.llh.T

    # Equation 7.19a and 7.19b from IERS Conventions 2010
    de = (
        coeff["A_d1_e"](lon, lat, grid=False) * np.cos(omega_1 * t)
        + coeff["B_d1_e"](lon, lat, grid=False) * np.sin(omega_1 * t)
        + coeff["A_d2_e"](lon, lat, grid=False) * np.cos(omega_2 * t)
        + coeff["B_d2_e"](lon, lat, grid=False) * np.sin(omega_2 * t)
    )
    dn = (
        coeff["A_d1_n"](lon, lat, grid=False) * np.cos(omega_1 * t)
        + coeff["B_d1_n"](lon, lat, grid=False) * np.sin(omega_1 * t)
        + coeff["A_d2_n"](lon, lat, grid=False) * np.cos(omega_2 * t)
        + coeff["B_d2_n"](lon, lat, grid=False) * np.sin(omega_2 * t)
    )
    du = (
        coeff["A_d1_u"](lon, lat, grid=False) * np.cos(omega_1 * t)
        + coeff["B_d1_u"](lon, lat, grid=False) * np.sin(omega_1 * t)
        + coeff["A_d2_u"](lon, lat, grid=False) * np.cos(omega_2 * t)
        + coeff["B_d2_u"](lon, lat, grid=False) * np.sin(omega_2 * t)
    )
    denu = np.vstack([de, dn, du]).T * Unit.mm2m

    if position.is_position(dset.site_pos):
        pos_correction = position.PositionDelta(denu, system="enu", ref_pos=dset.site_pos, time=dset.time)
    elif position.is_posvel(dset.site_pos):
        # set velocity to zero
        denu = np.concatenate((denu, np.zeros(denu.shape)), axis=1)
        pos_correction = position.PosVelDelta(denu, system="enu", ref_pos=dset.site_pos, time=dset.time)
    else:
        log.fatal(f"dset.site_pos{dset.default_field_suffix} is not a PositionArray or PosVelArray.")

    # Add center of mass corrections
    if use_cmc:
        # Equation (7.20) in [1]
        coeff_cmc = apriori.get("atmospheric_tides_cmc")
        cmc = (
            coeff_cmc["A1"][None, :] * np.cos(omega_1 * t)[:, None]
            + coeff_cmc["B1"][None, :] * np.sin(omega_1 * t)[:, None]
            + coeff_cmc["A2"][None, :] * np.cos(omega_2 * t)[:, None]
            + coeff_cmc["B2"][None, :] * np.sin(omega_2 * t)[:, None]
        )
        if position.is_position(dset.site_pos):
            cmc_correction = position.PositionDelta(cmc, system="trs", ref_pos=dset.site_pos, time=dset.time)
        elif position.is_posvel(dset.site_pos):
            # set velocity to zero
            cmc = np.concatenate((cmc, np.zeros(cmc.shape)), axis=1)
            cmc_correction = position.PosVelDelta(cmc, system="trs", ref_pos=dset.site_pos, time=dset.time)
        pos_correction = pos_correction.trs + cmc_correction.trs

    return pos_correction.gcrs
