"""Get apriori data for ephemeris

Description:

Reads data from CPF file



"""
import numpy as np

# Where imports
from where import parsers
from where.lib import config
from where.lib import files
from where.lib import log
from where.lib import plugins
from where.lib import time
from where.apriori import get_satellite_vars
from midgard.math import interpolation


@plugins.register
def get_ephemeris(rundate, sat_name):
    """Get CPF data for a given date

    The ephemeris data is stored in a dictionary with keys site_num (5 digits, e.g. '97401') and subdictionary with
    keys antenna_num (one letter and three digits, e.g. 'S002'). Contains information on site position and sigmas.

    Args:
        rundate (Datetime):   Model run date.
        sat_name (String):    Name of satellite.

    Returns:
        Dict: Ephemeris data.
    """
    file_key = "slr_ephemeris"
    sat_data = get_satellite_vars(sat_name)

    provider_list = config.tech.prediction_providers.list
    # Find the latest version of the observation file
    versions = files.glob_variable(file_key, "version", r"\d{4}", file_vars=sat_data)
    ephemeris_data = dict()

    try:
        ephemeris_data["version"] = sorted(versions)[-1]
        providers = files.glob_variable(file_key, "provider", r"\w+", file_vars=ephemeris_data)
        for provider in provider_list:
            if provider in providers:
                ephemeris_data["provider"] = provider
                break
    except IndexError:
        print(f"Pattern: '{files.path(file_key)}'")  # TODO: Because of format log does not print this properly
        log.fatal(f"No ephemeris data found")

    eph = parsers.parse_key(file_key, file_vars=ephemeris_data, rundate=rundate)
    eph = calculate_initial_values(eph)

    return eph


def calculate_initial_values(eph):
    """Computing initial values for position and velocity in GCRS system

    This is for later use in orbit integration, from tables in the prediction files.  Use a lagrange polynomial in
    order to interpolate in the tables.

    Args:
        eph:  Dict containing ephemeris information

    Returns:
        eph:  Dict where the initial position and velocity is added
    """
    pos_gcrs = np.empty((3, 0))
    times = np.empty((0))
    table_of_positions = sorted(eph.data["positions"].items())
    mjd1, mjd2 = zip(*[t for t, d in table_of_positions])

    for pos_time, (_, data) in zip(time.Time(val=mjd1, val2=mjd2, format="mjd", scale="utc"), table_of_positions):
        diffsec = (pos_time.utc.datetime - eph.rundate).total_seconds()
        # Only look at points close to rundate (start of integration)
        # if abs(diffsec) > 4000:
        #    continue
        # Table given in ITRF coordinate system. Convert to GCRS, where the integration of the satellite orbit will
        # be done
        pos_gcrs = np.hstack((pos_gcrs, np.transpose([pos_time.itrs2gcrs @ data["pos"]])))
        times = np.hstack((times, diffsec))

    log.info("Interpolating data from prediction file in order to get initial pos/vel")
    pos_gcrs_ip, vel_gcrs_ip = interpolation.interpolate_with_derivative(
        times, np.transpose(pos_gcrs), np.array([0.0]), kind="lagrange", window=10, bounds_error=False
    )
    eph["initial_pos"] = pos_gcrs_ip[0]
    eph["initial_vel"] = vel_gcrs_ip[0]
    return eph
