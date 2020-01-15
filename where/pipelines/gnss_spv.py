"""a GNSS pipeline

Description:
------------

TODO

"""

# Standard library imports
import itertools

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import gnss
from where.lib import log
from where import apriori
from where import cleaners
from where import estimation
from where.models import site, delay
from where import parsers
from where import writers


# The name of this technique
TECH = __name__.split(".")[-1]

TEST = False


@plugins.register_named("options")
def options():
    """Command line options that can be used to specify this technique

    Returns:
        Tuple:  Strings specifying command line options.
    """
    return ("--gnss_spv",)


@plugins.register_named("list_sessions")
def list_sessions(rundate):
    """Get list of sessions for the given rundate

    In the GNSS analysis the sessions correspond to stations. The list of stations is looked up in the config file.

    Args:
        rundate (date):   The model run date.

    Returns:
        List: Strings with names of stations.
    """
    return config.where[TECH].stations.list


@plugins.register_named("file_vars")
def file_vars():
    """File variables that will be available during the running of this technique

    In addition, date and analysis variables are available.

    Returns:
        Dict:  File variables special for this technique.
    """
    if not config.analysis.session.str:
        log.fatal(
            f"Session name is not defined. Use 'session' option for defining which station should be used in {TECH.upper()} analysis."
        )
    session = config.analysis.session.str
    return dict(station=session, STATION=session.upper())


#
# READ DATA
#
@plugins.register
def read(stage, dset):
    """Read the GNSS RINEX data.

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    dset.vars.update(file_vars())
    station = dset.vars["station"]

    # Read GNSS observation data either from Android raw file or RINEX file
    if config.tech.format.str == "android":
        parser = parsers.parse("gnss_android_raw_data", rundate=dset.analysis["rundate"], station=station)
    else:
        version, filepath = gnss.get_rinex_file_version("gnss_rinex_obs")
        log.info(f"RINEX file format {version}")
        if version.startswith("2"):
            parser = parsers.parse("rinex2_obs", rundate=dset.analysis["rundate"], station=station)
        elif version.startswith("3"):
            parser = parsers.parse("rinex3_obs", rundate=dset.analysis["rundate"], station=station)
        else:
            log.fatal(f"Unknown RINEX format {version} is used in file {filepath}")

    parser.write_to_dataset(dset)
    dset.write_as(stage=stage)
    dset.read()  # TODO: workaround because caching does not work correctly


#
# ORBIT DETERMINATION
#
@plugins.register
def orbit(stage, dset):
    """Determine GNSS satellite orbit

    TODO: Is the workflow for determining the satellite transmission time correct? gLAB determines satellite clock
          correction based on receiver time and not an satellite transmission time. Additionally gLAB does not apply
          relativistic corrections.

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    station = dset.vars["station"]
    orb_flag = config.tech.apriori_orbit.str

    # First estimate of satellite transmission time
    sat_time = dset.time - gnss.get_initial_flight_time(dset)

    # Second estimate using satellite clock and relativistic clock corrections
    orbit = apriori.get(
        "orbit", rundate=dset.analysis["rundate"], system=tuple(dset.unique("system")), station=station
    )
    orbit.dset_raw.write_as(stage=stage, session=station, dataset_name="raw")
    orbit.dset_edit.write_as(stage=stage, session=station, dataset_name="edit")

    # Determine initial satellite orbit solution with observation time as approximation
    orbit.calculate_orbit(dset, time="time")

    # TODO: Has it an effect to iterate here to improve satellite transmission time?
    # Determine satellite transmission time based on initial satellite orbit solution
    dset.add_time(
        "sat_time",
        val=dset.time
        - gnss.get_initial_flight_time(
            dset, sat_clock_corr=orbit.dset.gnss_satellite_clock, rel_clock_corr=orbit.dset.gnss_relativistic_clock
        ),
        scale=dset.time.scale,
    )

    # Use satellite transmission time for determination of satellite orbits
    orbit.calculate_orbit(dset, time="sat_time")

    # Copy to regular dataset
    dset.add_posvel("sat_posvel", itrs=orbit.dset.sat_posvel.itrs, time="sat_time", other="site_pos")
    # TODO: Is is possible to set the "calc_models" table in the model part?
    dset.add_float("gnss_satellite_clock", val=-orbit.dset.gnss_satellite_clock, table="calc_models")
    dset.add_float("gnss_relativistic_clock", val=-orbit.dset.gnss_relativistic_clock, table="calc_models")
    dset.site_pos.connect(dset.sat_posvel)
    dset.vars["orbit"] = orb_flag  # Needed e.g. for calling gnss_relativistic_clock model correctly.

    # Connect site position with satellite orbits needed for determination of elevation and azimuth
    dset.site_pos.connect(dset.sat_posvel)

    # Correct satellite position/velocity due Earth's rotation effect during signal flight time
    sat_pos, sat_vel = gnss.get_earth_rotation(dset)
    log.warn(
        "Correction of satellite position/velocity due to Earth's rotation effect during signal flight time is not"
        " applied."
    )
    # dset.sat_posvel.add_to_itrs(
    #    np.hstack((sat_pos, sat_vel))
    # )  # TODO: Check residuals will be worser by using that. Why?
    dset.add_posvel("gnss_earth_rotation", itrs=np.hstack((sat_pos, sat_vel)), time="sat_time")
    dset.write_as(stage=stage)
    dset.read()  # TODO: workaround because caching does not work correctly


#
# EDIT DATA
#
@plugins.register
def edit(stage, dset):
    """Edit GNSS data

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    cleaners.apply_editors("editors", dset)
    cleaners.apply_removers("removers", dset)
    dset.write_as(stage=stage)
    dset.read()  # TODO: workaround because caching does not work correctly


#
# CALCULATE AND ESTIMATE
#
# @plugins.register
def calculate_estimate(stage, dset):
    """Calculate model parameters and estimate

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    max_iterations = config.tech.max_iterations.int

    for iter_num in itertools.count(start=1):

        # CALCULATE
        # -----------
        # Correction of station position in GCRS due to loading and tide effects
        site.calculate_site("site", dset, shape=(3,))
        delta_pos = np.sum(dset.get_table("site").reshape((dset.num_obs, -1, 3)), axis=1)
        dset.site_pos.add_to_gcrs(delta_pos)

        # Initialize models given in configuration file by adding model fields to Dataset
        delay.calculate_delay("calc_models", dset, write_levels=dict(gnss_range="operational"))
        if "obs" in dset.fields:
            dset.obs[:] = gnss.get_code_observation(dset)
        else:
            dset.add_float("obs", val=gnss.get_code_observation(dset), unit="meter")

        # Get model corrections
        if "calc" in dset.fields:
            dset.calc[:] = np.sum(dset.get_table("calc_models"), axis=1)
        else:
            dset.add_float("calc", val=np.sum(dset.get_table("calc_models"), axis=1), unit="meter")

        if "residual" in dset.fields:
            dset.residual[:] = dset.obs - dset.calc
        else:
            dset.add_float("residual", val=dset.obs - dset.calc, unit="meter")

        # Store calculate results
        log.info(f"{dset.num_obs} observations, residual = {dset.rms('residual'):.4f}")
        dset.write_as(stage="calculate", dataset_id=iter_num)
        dset.read()  # TODO: workaround because caching does not work correctly

        # ESTIMATE
        # ----------
        partial_vectors = estimation.partial_vectors(dset, "estimate_method")

        log.blank()  # Space between iterations for clarity
        log.info(f"Estimating parameters for iteration {iter_num}")
        estimation.call("estimate_method", dset=dset, partial_vectors=partial_vectors, obs_noise=np.ones(dset.num_obs))
        rms = dset.rms("residual")
        log.info(f"{dset.num_obs} observations, postfit residual = {rms:.4f}")

        dset.write_as(stage="estimate", dataset_id=iter_num - 1)
        dset.read()  # TODO: workaround because caching does not work correctly

        # Detect and remove outliers based on residuals
        keep_idx = estimation.detect_outliers("estimate_outlier_detection", dset)

        if dset.meta["estimate_convergence_status"] and keep_idx.all():
            log.info(f"Estimation convergence limit of {config.tech.convergence_limit.float:.3e} is fulfilled.")
            break
        if iter_num >= max_iterations:
            break

        dset.subset(keep_idx)
        log.blank()


# ===========================================
# 	Single Point Velocity by Doppler
#   		(ECEF = ITRS)
# ===========================================
@plugins.register
def spv_doppler(stage, dset):
    """Calculate model parameters and estimate

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """

    # to have access to orbit data, you have make a new call
    # to orbit
    # orbit = apriori.get("orbit", rundate=dset.analysis["rundate"], system=tuple(dset.unique('system')), station=dset.vars["station"],)

    station = dset.vars["station"]
    orbit = apriori.get(
        "orbit", rundate=dset.analysis["rundate"], system=tuple(dset.unique("system")), station=station
    )
    orbit.dset_raw.write_as(stage=stage, session=station, dataset_name="raw")
    orbit.dset_edit.write_as(stage=stage, session=station, dataset_name="edit")

    import IPython

    IPython.embed()

    #  get the observed time
    MJD = sorted(set(dset.time.gps.mjd))
    for epoch in MJD:
        idx = dset.time.gps.mjd == epoch
        valid_sats = dset.satellite[idx]
        dop_obs = dset.D1X[idx]
        AZ = dset.sat_posvel.azimuth[idx]
        EL = dset.sat_posvel.elevation[idx]
        SAT_ENU = dset.sat_posvel.enu_up[idx]
        SAT_POS = dset.sat_posvel.itrs_pos[idx]
        SAT_VEL = dset.sat_posvel.itrs_vel[idx]

        # valid_sats = pd.DataFrame(dset.satellite[idx])
        # dop_obs    = pd.DataFrame(dset.D1C[idx]
        # my_DF      = pd.concat([valid_sats, dop_obs], axis=1)


#
# WRITE RESULTS
#
@plugins.register
def write(stage, dset):
    """Write results to file.

    Write results to file. This uses the writers framework which calls different writers depending on the output-field
    in the config-file.

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    writers.write(default_dset=dset)
