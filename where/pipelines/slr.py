"""a SLR pipeline

Description:
------------

Calls the various models and perform the necessary SLR computations.

"""
# Standard library imports
from datetime import datetime, timedelta
import itertools
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math import interpolation
from midgard.math.constant import constant

# Where imports
from where import apriori
from where import cleaners
from where import estimation
from where.models import site, delay, orbit
from where import obs
from where import writers
from where.lib import config
from where.lib import log
from where.lib import rotation
from where.data.time import Time, TimeDelta


# The name of this technique
TECH = __name__.split(".")[-1]


@plugins.register_named("options")
def options():
    """Command line options that can be used to specify this technique

    Returns:
        Tuple:  Strings specifying command line options.
    """
    return "-s", "--slr"


@plugins.register_named("list_sessions")
def list_sessions():
    """Get list of sessions for the given rundate

    In the SLR analysis the sessions correspond to satellites. The list of satellites is looked up in the config file.

    Returns:
        List: Strings with names of satellites.
    """
    #    return config.tech.satellites.list
    return [""]


@plugins.register_named("validate_args")
def validate_args(rundate, **kwargs):
    """Validate a session

    No validation needed.

    Args:
        rundate:           the date we run the analysis for
        kwargs:            keyword arguments from command line

    """


@plugins.register_named("file_vars")
def file_vars(file_vars=None):
    """File variables that will be available during the running of this technique

    In addition, date and analysis variables are available.

    Returns:
        Dict:  File variables special for this technique.
    """
    _file_vars = dict()

    # Sinex file vars
    if "sinex" in config.tech.section_names:
        _file_vars["solution"] = config.tech.sinex.solution.str
        _file_vars["file_agency"] = config.tech.sinex.file_agency.str.lower()

    return _file_vars


@plugins.register
def read(stage, dset):
    """Read the SLR data

    Args:
        stage:  Name of current stage
        dset:   Dataset containing the data
    """
    sat_name = dset.vars["sat_name"]
    obs.get(dset, sat_name=sat_name)
    log.info(f"Parsed {dset.num_obs} observations")
    dset.write_as(stage=stage, label=0)


@plugins.register
def edit(stage, dset):
    """Edit the data by applying editor

    Args:
        stage:  Name of current stage
        dset:   Dataset containing the data
    """
    # Clean up dataset
    cleaners.apply_removers("removers", dset)

    # Indicate if range and time bias are estimated or not
    # and apply biases
    dset.add_float("range_bias", np.zeros(dset.num_obs), unit="meter")
    dset.add_float("time_bias", np.zeros(dset.num_obs), unit="meter")
    dset.add_bool("estimate_range", np.zeros(dset.num_obs))
    dset.add_bool("estimate_time", np.zeros(dset.num_obs))

    for station in config.tech.slr_range_bias.estimate_stations.list:
        int_idx = dset.filter(station=station)
        if np.any(int_idx):
            log.info(f"Config file: Will estimate range bias for station {station} in estimation stage")
            dset.estimate_range[:] = np.logical_or(int_idx, dset.estimate_range[:])
    for station in config.tech.slr_time_bias.estimate_stations.list:
        int_idx = dset.filter(station=station)
        if np.any(int_idx):
            log.info(f"Config file: Will estimate time bias for station {station} in estimation stage")
            dset.estimate_time[:] = np.logical_or(int_idx, dset.estimate_time[:])

    cleaners.apply_editors("editors", dset)

    # Write dataset
    dset.write_as(stage=stage, label=0)


@plugins.register
def calculate(stage, dset):
    """
    Integrate differential equation of motion of the satellite

    Args:
        stage:  Name of current stage
        dset:   Dataset containing the data
    """

    iterations = config.tech.iterations.int

    # Run models adjusting station positions
    site.calculate_site("site", dset)
    delta_pos = site.add("site", dset)
    dset.site_pos[:] = (dset.site_pos.gcrs + delta_pos[0].gcrs).trs

    dset.add_float("obs", val=dset.time_of_flight * constant.c / 2, unit="meter")
    dset.add_float("calc", np.zeros(dset.num_obs), unit="meter")
    dset.add_float("residual", np.zeros(dset.num_obs), unit="meter")
    dset.add_float("up_leg", np.zeros(dset.num_obs), unit="second")
    dset.add_posvel("sat_pos", np.zeros((dset.num_obs, 6)), system="gcrs", time=dset.time)
    arc_length = config.tech.arc_length.float

    dset.site_pos.other = dset.sat_pos

    # First guess for up_leg:
    dset.up_leg[:] = dset.time_of_flight / 2

    for iter_num in itertools.count(start=1):
        log.blank()
        log.info(f"Calculating model corrections for iteration {iter_num}")

        sat_time_list = dset.obs_time + dset.time_bias + dset.up_leg
        apriori_orbit_provider = config.tech.apriori_orbit.str
        sat_name = dset.vars["sat_name"]

        rundate = dset.analysis["rundate"]

        if apriori_orbit_provider:
            version = config.tech.apriori_orbit_version.str
            log.info(f"Using external orbits from {apriori_orbit_provider}, version {version}")
            apriori_orbit = apriori.get(
                "orbit",
                rundate=rundate + timedelta(days=arc_length),
                time=None,
                day_offset=6,
                satellite=sat_name,
                apriori_orbit="slr",
                file_key="slr_external_orbits",
            )
            dset_external = apriori_orbit._read(dset, apriori_orbit_provider, version)

            sat_pos = dset_external.sat_pos.gcrs_pos
            t_sec = TimeDelta(
                dset_external.time
                - Time(datetime(rundate.year, rundate.month, rundate.day), scale="utc", fmt="datetime"),
                fmt="seconds",
            )
            t_sec = t_sec.value
        else:
            sat_pos, sat_vel, t_sec = orbit.calculate_orbit(
                datetime(rundate.year, rundate.month, rundate.day), sat_name, sat_time_list, return_full_table=True
            )

        sat_pos_ip, sat_vel_ip = interpolation.interpolate_with_derivative(
            np.array(t_sec), sat_pos, sat_time_list, kind="interpolated_univariate_spline"
        )
        dset.sat_pos.gcrs[:] = np.concatenate((sat_pos_ip, sat_vel_ip), axis=1)
        delay.calculate_delay("kinematic_models", dset)

        # We observe the time when an observation is done, and the time of flight of the laser pulse. We estimate
        # the up-leg time with Newton's method applied to the equation (8.84) of :cite:'beutler2005' Gerhard Beutler:
        # Methods of Celestial Mechanics, Vol I., 2005.
        for j in range(0, 4):
            reflect_time = dset.time + TimeDelta(dset.time_bias + dset.up_leg, fmt="seconds", scale="utc")
            site_pos_reflect_time = (rotation.trs2gcrs(reflect_time) @ dset.site_pos.trs.val[:, :, None])[:, :, 0]
            sta_sat_vector = dset.sat_pos.gcrs.pos.val - site_pos_reflect_time
            unit_vector = sta_sat_vector / np.linalg.norm(sta_sat_vector, axis=1)[:, None]

            rho12 = (np.linalg.norm(sta_sat_vector, axis=1) + delay.add("kinematic_models", dset)) / constant.c
            correction = (-dset.up_leg + rho12) / (
                np.ones(dset.num_obs) - np.sum(unit_vector / constant.c * dset.sat_pos.vel.val, axis=1)
            )
            dset.up_leg[:] += correction
            sat_time_list = dset.obs_time + dset.time_bias + dset.up_leg
            sat_pos_ip, sat_vel_ip = interpolation.interpolate_with_derivative(
                np.array(t_sec), sat_pos, sat_time_list, kind="interpolated_univariate_spline"
            )

            dset.sat_pos.gcrs[:] = np.concatenate((sat_pos_ip, sat_vel_ip), axis=1)

        delay.calculate_delay("satellite_models", dset)
        dset.calc[:] = delay.add("satellite_models", dset)
        dset.residual[:] = dset.obs - dset.calc
        log.info(f"{dset.num_obs} observations, residual = {dset.rms('residual'):.4f}")
        if not apriori_orbit_provider:
            orbit.update_orbit(
                sat_name, dset.site_pos.gcrs, dset.sat_pos.pos, dset.sat_pos.vel, dset.residual, dset.bin_rms
            )

        dset.write_as(stage=stage, label=iter_num, sat_name=sat_name)
        if iter_num >= iterations:
            break


@plugins.register
def estimate(stage, dset):
    """Filter residuals

    Args:
        prev_stage (String): Name of previous stage.
        stage (String):      Name of current stage.
    """

    partial_vectors = estimation.partial_vectors(dset, "estimate_method")
    max_iterations = config.tech.estimate_max_iterations.int

    for iter_num in itertools.count(start=1):
        log.info(f"Estimating parameters for iteration {iter_num}")
        estimation.call(
            "estimate_method", dset=dset, partial_vectors=partial_vectors, obs_noise=None  # TODO: Add something here
        )
        rms = dset.rms("residual")
        log.info(f"{dset.num_obs} observations, postfit residual = {rms:.4f}")
        dset.write_as(stage=stage, label=iter_num - 1, sat_name=dset.vars["sat_name"])

        # TODO:
        # Do some iteration with removal of outliers?
        break

        if iter_num >= max_iterations:
            break

    estimation.solve_neq(dset)
    dset.write()


@plugins.register
def write(stage, dset):
    """Write results to file

    Write results to file. This uses the writers framework which calls different writers depending on the output-field
    in the config-file.

    Args:
        stage (String):      Name of current stage
        dset:                Dataset where data is stored.
    """

    writers.write(dset)
