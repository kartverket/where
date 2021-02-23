"""Get apriori data for ephemeris

Description:

Reads data from CPF file



"""
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math import interpolation

# Where imports
from where import parsers
from where.lib import config
from where.lib import log
from where.lib import rotation
from where.data import time
from where.apriori import get_satellite_vars


@plugins.register
def get_ephemeris(rundate, sat_name):
    """Get CPF data for a given date

    The ephemeris data is stored in a dictionary with tables of times and position under the "positions" key.
    This table is interpolated in the calculate_initial_values method.
    Args:
        rundate (Datetime):   Model run date.
        sat_name (String):    Name of satellite.

    Returns:
        Dict: Ephemeris data.
    """
    file_key = "slr_ephemeris"
    ephemeris_data = get_satellite_vars(sat_name)
    provider_list = config.tech.prediction_providers.list
    # Find the latest version of the observation file
    versions = config.files.glob_variable(file_key, "version", r"\d{4}", file_vars=ephemeris_data)

    try:
        ephemeris_data["version"] = sorted(versions)[-1]
        providers = config.files.glob_variable(file_key, "provider", r"\w+", file_vars=ephemeris_data)
        for provider in provider_list:
            if provider in providers:
                ephemeris_data["provider"] = provider
                break
        else:
            log.fatal(f"No valid provider found: {', '.join(providers)}")
    except IndexError:
        log.info("No ephemeris data found")
        log.info(f"Download manually from https://cddis.nasa.gov/archive/slr/cpf_predicts/{rundate.year}/{sat_name}")
        log.fatal(f"Please save missing file as '{config.files.path(file_key)}' !")
    eph_parser = parsers.parse_key(file_key, file_vars=ephemeris_data)
    eph = calculate_initial_values(eph_parser.as_dict(), rundate)

    return eph


def calculate_initial_values(eph, rundate):
    """Computing initial values for position and velocity in GCRS system

    This is for later use in orbit integration, from tables in the prediction files.  Use a lagrange polynomial in
    order to interpolate in the tables.

    Args:
        eph:  Dict containing ephemeris information

    Returns:
        eph:  Dict where the initial position and velocity is added
    """
    data = sorted(eph["positions"].items())
    pos_itrs = np.zeros((len(data), 3))
    mjd1, mjd2 = zip(*[t for t, d in data])
    rotation_mat = rotation.trs2gcrs(time.Time(val=mjd1, val2=mjd2, fmt="mjd", scale="utc"))
    tbl = time.Time(val=mjd1, val2=mjd2, fmt="mjd", scale="utc")

    for i in range(0, len(data)):
        pos_itrs[i] = data[i][1]["pos"]

    diffsec = np.array([(t - rundate).total_seconds() for t in tbl.utc.datetime])

    # Table given in ITRF coordinate system. Convert to GCRS, where the integration of the satellite orbit will
    # be done

    pos_gcrs = np.sum(rotation_mat @ pos_itrs[:, :, None], axis=2)
    log.info("Interpolating data from prediction file in order to get initial pos/vel")
    pos_gcrs_ip, vel_gcrs_ip = interpolation.interpolate_with_derivative(
        diffsec, pos_gcrs, np.array([0.0]), kind="lagrange", window=10, bounds_error=False
    )
    eph["initial_pos"] = pos_gcrs_ip[0]
    eph["initial_vel"] = vel_gcrs_ip[0]

    return eph
