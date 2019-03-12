"""a VLBI pipeline

Description:
------------

TODO

"""
import itertools

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where import obs
from where import cleaners
from where import estimation
from where.lib import config
from where.lib import exceptions
from where.lib import files
from where.lib import log
from where.lib import util
from where import models
from where import writers

# The name of this technique
TECH = __name__.split(".")[-1]


@plugins.register_named("options")
def options():
    """Command line options that can be used to specify this technique

    Returns:
        Tuple:  Strings specifying command line options.
    """
    return "-v", "--vlbi"


@plugins.register_named("list_sessions")
def list_sessions(rundate):
    """Sessions available for the given rundate

    Args:
        rundate (date):   The model run date.

    Returns:
        List:   Strings with names of available sessions.
    """
    if config.where.get(
        "get_session_from_master",
        section=TECH,
        value=util.read_option_value("--get_session_from_master", default=None),  # TODO: add this to mg_config
        default=False,
    ).bool:
        skip_sessions = set(
            config.where.get(
                "skip_sessions",
                section="runner",
                value=util.read_option_value("--skip_sessions", default=None),
                default="",
            ).list
        )
        session_types = config.where.get(
            "session_types",
            section="runner",
            value=util.read_option_value("--session_types", default=None),
            default="",
        ).list
        master_schedule = apriori.get("vlbi_master_schedule", rundate=rundate)
        sessions = set(master_schedule.list_sessions(rundate, session_types=session_types))
        sessions = sessions - skip_sessions

        return sessions
    else:
        obs_format = config.tech.get("obs_format", section=TECH).str  # TODO: This always falls back on config.where ..
        file_vars = config.create_file_vars(rundate, TECH, session=None)
        del file_vars["session"]  # TODO: Do not add None variables to file_vars?
        found_sessions = files.glob_variable(
            f"vlbi_obs_{obs_format}", variable="session", pattern=r"\w{2}", file_vars=file_vars
        )
        return found_sessions


@plugins.register_named("validate_session")
def validate_session(rundate, session):
    """Validate a session for the given rundate

    If session is not a valid VLBI session for the given rundate, an InvalidSessionError is raised.

    Args:
        rundate (date):    The model run date.
        session (String):  Name of session.

    Return:
        String:  Name of validated session.
    """
    if not session:
        if util.check_options("--session"):
            return session  # Explicitly specified blank session, typically to open timeseries interactively
        raise exceptions.InvalidSessionError("You must specify '--session=<...>' to run a VLBI analysis")

    # TODO: Can we use master files to validate sessions? What about intensives?
    master_schedule = apriori.get("vlbi_master_schedule", rundate=rundate)
    master_sessions = master_schedule.list_sessions(rundate)
    if session not in master_sessions:
        log.warn(
            f"Session '{session}' is not listed in master file for {rundate:{config.FMT_date}}. "
            f"Available sessions are {', '.join(master_sessions)}."
        )
    return session


@plugins.register_named("file_vars")
def file_vars():
    """File variables that will be available during the running of this technique

    In addition, date and analysis variables are available.

    Returns:
        Dict:  File variables special for this technique.
    """
    file_vars = dict()

    # Add obs_version for ngs
    if config.tech.get("obs_format").str == "ngs":
        versions = files.glob_variable("vlbi_obs_ngs", "obs_version", r"\d{3}")
        if versions:
            file_vars["obs_version"] = max(versions)
        elif config.where.files.download_missing.bool:
            # Look online for a candidate
            log.info("No NGS observation file found on disk: Looking for one online.")
            obs_versions = [f"{v:03d}" for v in reversed(range(4, 10))]
            for obs_version in obs_versions:
                url = files.url(
                    "vlbi_obs_ngs", file_vars=dict(obs_version=obs_version), is_zipped=True, use_aliases=False
                )
                log.info(f"Looking for {url} ...")
                if url.exists():
                    file_vars["obs_version"] = obs_version
                    break

        if not file_vars:
            log.fatal("No NGS observation file found")

    # Add obs_version for vgosdb
    if config.tech.get("obs_format").str == "vgosdb":
        versions = files.glob_variable("vlbi_obs_vgosdb", "obs_version", r"\d{3}")
        if versions:
            file_vars["obs_version"] = max(versions)
        elif config.where.files.download_missing.bool:
            # Look online for a candidate
            log.warn("No VGOSDB wrapper file found. Not attempting to download. TODO")
            # log.info("No NGS wrapper file found on disk: Looking for one online.")
            # obs_versions = [f"{v:03d}" for v in reversed(range(4, 10))]
            # for obs_version in obs_versions:
            #    url = files.url(
            #        "vlbi_obs_ngs", file_vars=dict(obs_version=obs_version), is_zipped=True, use_aliases=False
            #    )
            #    log.info(f"Looking for {url} ...")
            #    if url.exists():
            #        file_vars["obs_version"] = obs_version
            #        break

        if not file_vars:
            log.fatal("No VGOSDB observation file found")

    # Sinex file vars
    if "sinex" in config.tech.section_names:
        file_vars["solution"] = config.tech.sinex.solution.str
        file_vars["file_agency"] = config.tech.sinex.file_agency.str.lower()

    return file_vars


@plugins.register
def read(stage, dset):
    """Read VLBI data

    Args:
        dset (Dataset):  Dataset
    """
    obs.get(dset)
    log.info(f"Parsed {dset.num_obs} observations")

    dset.write_as(stage=stage, dataset_id=0)
    dset.read()


@plugins.register
def edit(stage, dset):
    """Edit the observations

    Args:
        rundate (Datetime):  The model run date.
        session (String):    Name of session.
        prev_stage (String): Name of previous stage.
        stage (String):      Name of current stage.
    """
    cleaners.apply_editors("editors", dset)
    cleaners.apply_removers("removers", dset)
    dset.write_as(stage=stage, dataset_id=0)
    dset.read()


@plugins.register
def calculate(stage, dset):
    """Estimate model parameters

    Args:
        rundate (Datetime):  The model run date.
        session (String):    Name of session.
        prev_stage (String): Name of previous stage.
        stage (String):      Name of current stage.
    """
    # Run models adjusting station positions
    log.info(f"Calculating station displacements")
    models.calculate_site("pos_models", dset, shape=(6,))
    delta_pos = np.sum(dset.get_table("pos_models").reshape((dset.num_obs, -1, 6)), axis=1)
    gcrs_dpos_1 = delta_pos[:, :3]
    gcrs_dvel_1 = (dset.time.itrs2gcrs_dot @ dset.site_pos_1.convert_gcrs_to_itrs(gcrs_dpos_1)[:, :, None])[:, :, 0]
    dset.site_pos_1.add_to_gcrs(np.concatenate((gcrs_dpos_1, gcrs_dvel_1), axis=1))
    gcrs_dpos_2 = delta_pos[:, 3:]
    gcrs_dvel_2 = (dset.time.itrs2gcrs_dot @ dset.site_pos_2.convert_gcrs_to_itrs(gcrs_dpos_2)[:, :, None])[:, :, 0]
    dset.site_pos_2.add_to_gcrs(np.concatenate((gcrs_dpos_2, gcrs_dvel_2), axis=1))
    log.blank()

    # Run models for each term of the observation equation
    log.info(f"Calculating theoretical delays")
    models.calculate_delay("calc_models", dset)
    dset.add_float("obs", val=dset.observed_delay, unit="meter", write_level="operational")
    dset.add_float("calc", val=np.sum(dset.get_table("calc_models"), axis=1), unit="meter", write_level="operational")
    dset.add_float("residual", val=dset.obs - dset.calc, unit="meter", write_level="operational")
    log.blank()

    # Estimate clock polynomial
    log.info(f"Calculating clock polynomials")
    max_iterations = config.tech.calculate_max_iterations.int
    outlier_limit = config.tech.calculate_outlier_limit.float
    store_outliers = config.tech.store_outliers.bool

    for iter_num in itertools.count(start=1):
        models.calculate_delay("correction_models", dset, dset)
        dset.calc[:] = np.sum(np.hstack((dset.get_table("calc_models"), dset.get_table("correction_models"))), axis=1)
        dset.residual[:] = dset.obs - dset.calc
        rms = dset.rms("residual")
        log.info(f"{dset.num_obs} observations, residual = {rms:.4f}")

        # Store results
        dset.write_as(stage=stage, dataset_id=iter_num - 1)
        dset.read()

        # Detect and remove extreme outliers
        idx = np.abs(dset.residual) < outlier_limit * rms
        if iter_num > max_iterations or idx.all():
            break

        if store_outliers:
            bad_idx = np.logical_not(idx)
            log.info(f"Adding {np.sum(bad_idx)} observations to ignore_observation")
            bad_obs = np.char.add(np.char.add(dset.time.utc.iso[bad_idx], " "), dset.baseline[bad_idx]).tolist()
            with config.update_tech_config(dset.rundate, TECH, dset.vars["session"]) as cfg:
                current = cfg.ignore_observation.observations.as_list(", *")
                updated = ", ".join(sorted(current + bad_obs))
                cfg.update("ignore_observation", "observations", updated, source=util.get_program_name())

        dset.subset(idx)
        log.info(f"Removing {sum(~idx)} observations with residuals bigger than {outlier_limit * rms}")
        log.blank()

    # Try to detect clock breaks
    if config.tech.detect_clockbreaks.bool:
        writers.write_one("vlbi_detect_clockbreaks", dset=dset)
        dset.write()


@plugins.register
def estimate(stage, dset):
    """Filter residuals

    Args:
        rundate (Datetime):  The model run date.
        session (String):    Name of session.
        prev_stage (String): Name of previous stage.
        stage (String):      Name of current stage.
    """
    partial_vectors = estimation.partial_vectors(dset, "estimate_method")
    max_iterations = config.tech.estimate_max_iterations.int

    for iter_num in itertools.count(start=1):
        log.info(f"Estimating parameters for iteration {iter_num}")
        estimation.call(
            "estimate_method",
            dset=dset,
            partial_vectors=partial_vectors,
            obs_noise=dset.observed_delay_ferr ** 2 + 0.01 ** 2,
        )
        rms = dset.rms("residual")
        log.info(f"{dset.num_obs} observations, postfit residual = {rms:.4f}")
        dset.write_as(stage=stage, dataset_id=iter_num - 1)
        dset.read()
        if iter_num >= max_iterations:
            break

        # Detect and remove outliers
        keep_idx = estimation.detect_outliers("estimate_outlier_detection", dset)
        if keep_idx.all():
            break
        dset.subset(keep_idx)
        log.blank()


@plugins.register
def write(stage, dset):
    """Write results to file

    Write results to file. This uses the writers framework which calls different writers depending on the output-field
    in the config-file.

    Args:
        rundate (Datetime):  The model run date.
        session (String):    Name of session.
        prev_stage (String): Name of previous stage.
        stage (String):      Name of current stage.
    """
    writers.write(default_dset=dset)
