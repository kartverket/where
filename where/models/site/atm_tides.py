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
from midgard.dev import exceptions
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where import apriori
from where.lib import config
from where.lib import log
from where.data import position


MODEL = __name__.split(".")[-1]

@plugins.register
def atm_tides(dset):
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
    use_cmc = config.tech[MODEL].cmc.bool

    # Time argument is fraction of UT1 day, see [2].
    t = dset.time.ut1.jd_frac

    coeff = get_coefficients(dset)

    # S1 has a period of 1 cycle/day, S2 has a period of 2 cycle/day
    omega_1 = 2 * np.pi
    omega_2 = 4 * np.pi

    # Equation 7.19a and 7.19b from IERS Conventions 2010
    de = (coeff["A_d1_e"] * np.cos(omega_1 * t)
        + coeff["B_d1_e"] * np.sin(omega_1 * t)
        + coeff["A_d2_e"] * np.cos(omega_2 * t)
        + coeff["B_d2_e"] * np.sin(omega_2 * t)
    )
    dn = (coeff["A_d1_n"] * np.cos(omega_1 * t)
        + coeff["B_d1_n"] * np.sin(omega_1 * t)
        + coeff["A_d2_n"] * np.cos(omega_2 * t)
        + coeff["B_d2_n"] * np.sin(omega_2 * t)
    )
    du = (coeff["A_d1_u"] * np.cos(omega_1 * t)
        + coeff["B_d1_u"] * np.sin(omega_1 * t)
        + coeff["A_d2_u"] * np.cos(omega_2 * t)
        + coeff["B_d2_u"] * np.sin(omega_2 * t)
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

def get_coefficients(dset):
    model_type = config.tech[MODEL].type.str

    stations = dset.unique("station")
    coefficients = [
        "A_d1_e", "B_d1_e", "A_d2_e", "B_d2_e",
        "A_d1_n", "B_d1_n", "A_d2_n", "B_d2_n",
        "A_d1_u", "B_d1_u", "A_d2_u", "B_d2_u",
        ]
    output = dict()
    for c in coefficients:
        output[c] = np.zeros(dset.num_obs)

    def _grid_coeff(stations, **kwargs):
        coeff_interpolator = apriori.get("atmospheric_tides_grid", **kwargs)

        for sta in stations:
            sta_idx = dset.filter(station=sta)
            lat, lon, _ = dset.site_pos[sta_idx][0].pos.llh
            for c in coefficients:
                output[c][sta_idx] = coeff_interpolator[c](lon, lat, grid=False)

    if model_type == "grid":
        _grid_coeff(stations)

    elif model_type == "station":
        model = config.tech[MODEL].model.str
        origin = config.tech[MODEL].origin.str
        missing_sta = set()
        try:
            file_key = f"atmospheric_tides_station_{model}_{origin}"
            coeff_sta = apriori.get(file_key)

            for sta in stations:
                sta_idx = dset.filter(station=sta)
                for c in coefficients:
                    try:
                        output[c][sta_idx] = coeff_sta[sta][c]
                    except KeyError:
                        missing_sta.add(sta)
        except exceptions.MissingSectionError:
            log.error(f"Unknown model. No file_key found for {file_key}")
        if missing_sta:
            log.warn(f"Missing APL coefficients for {', '.join(missing_sta)} for {file_key}. Using grid_{model}_ce as backup.")
            _grid_coeff(missing_sta, model=model, origin="ce")

    return output
