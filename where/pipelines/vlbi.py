"""a VLBI pipeline

Description:
------------


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
from where.lib import log
from where.lib import util
from where.models import site, delay
from where import writers

# The name of this technique
pipeline = __name__.split(".")[-1]


@plugins.register_named("options")
def options():
    """Command line options that can be used to specify this technique

    Returns:
        Tuple:  Strings specifying command line options.
    """
    return "-v", "--vlbi"


@plugins.register_named("get_args")
def get_args(rundate, input_args=None):
    """Convert where_runner arguments to where arguments for given date

    Args:
        rundate (date):   The model run date.

    Returns:
        List:   Strings with names of available sessions.
    """
    keyword = "--session"
    session_list = set()
    input_args = list(input_args) if input_args is not None else list()
    for idx in range(len(input_args)):
        key, _, value = input_args[idx].partition("=")
        if key == keyword:
            session_list = set(value.split(","))
            input_args.pop(idx)
            break
    args = " ".join(input_args)

    get_session_from_master = config.where.get(
        "get_session_from_master",
        section=pipeline,
        value=util.read_option_value("--get_session_from_master", default=None),  # TODO: add this to mg_config
        default=False,
    ).bool

    if get_session_from_master:
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

        check_master_status = config.where.get(
            "check_master_status",
            section="runner",
            value=util.read_option_value("--check_master_status", default=None),
            default=False,
        ).bool

        not_ready_sessions = set()
        if check_master_status:
            for session in sessions:
                if not master_schedule.ready(rundate, session):
                    status = master_schedule.status(rundate, session)
                    log.warn(
                        f"{rundate} {session} is not ready for processing. Master file status: '{status}'. Skipping session."
                    )
                    not_ready_sessions.add(session)

        sessions = set(sessions) - skip_sessions - not_ready_sessions
        sessions = sessions & session_list if session_list else sessions
        return [keyword + "=" + s + " " + args for s in sessions]
    else:
        obs_format = config.tech.get(
            "obs_format", section=pipeline
        ).str  # TODO: This always falls back on config.where ..
        file_vars = config.create_file_vars(rundate, pipeline, session=None)
        del file_vars["session"]  # TODO: Do not add None variables to file_vars?
        sessions = config.files.glob_variable(
            f"vlbi_obs_{obs_format}", variable="session", pattern=r"\w{2}", file_vars=file_vars
        )
        sessions = sessions & session_list
        return [keyword + "=" + s + " " + args for s in sessions]


@plugins.register_named("validate_args")
def validate_args(rundate, session=None, **kwargs):
    """Validate a session for the given rundate

    If session is not a valid VLBI session for the given rundate, an InvalidSessionError is raised.

    Args:
        rundate (date):    The model run date.
        session (String):  Name of session.

    Return:
        String:  Name of validated session.
    """
    if session is None:
        raise exceptions.InvalidArgsError("You must specify '--session=<...>' to run a VLBI analysis")

    master_schedule = apriori.get("vlbi_master_schedule", rundate=rundate)
    master_sessions = master_schedule.list_sessions(rundate)
    if session not in master_sessions:
        log.warn(
            f"Session '{session}' is not listed in master file for {rundate:{config.FMT_date}}. "
            f"Available sessions are {', '.join(master_sessions)}."
        )
    else:
        if not master_schedule.ready(rundate, session):
            status = master_schedule.status(rundate, session)
            log.warn(f"Session '{session}' is not ready for processing. Master file status: '{status}'")


@plugins.register_named("file_vars")
def file_vars():
    """File variables that will be available during the running of this technique

    In addition, date and analysis variables are available.

    Returns:
        Dict:  File variables special for this technique.
    """
    _file_vars = dict()

    # Add obs_version for ngs
    if config.tech.get("obs_format").str == "ngs":
        versions = config.files.glob_variable("vlbi_obs_ngs", "obs_version", r"\d{3}")
        if versions:
            _file_vars["obs_version"] = max(versions)
        elif config.where.files.download_missing.bool:
            # Look online for a candidate
            log.info("No NGS observation file found on disk: Looking for one online.")
            obs_versions = [f"{v:03d}" for v in reversed(range(4, 10))]
            for obs_version in obs_versions:
                url = config.files.url(
                    "vlbi_obs_ngs", file_vars=dict(obs_version=obs_version), is_zipped=True, use_aliases=False
                )
                log.info(f"Looking for {url} ...")
                if url.exists():
                    _file_vars["obs_version"] = obs_version
                    break

        if not _file_vars:
            log.fatal("No NGS observation file found")

    # Add obs_version for vgosdb
    if config.tech.get("obs_format").str == "vgosdb":
        versions = config.files.glob_variable("vlbi_obs_vgosdb", "obs_version", r"\d{3}")
        if versions:
            _file_vars["obs_version"] = max(versions)
        agencies = config.files.glob_variable("vlbi_obs_vgosdb", "agency", r"[\w]+", file_vars=_file_vars)
        if agencies:
            _file_vars["agency"] = "IVS" if "IVS" in agencies else agencies.pop()
            if len(agencies) > 1:
                log.warn(
                    f"Multiple agencies found ({', '.join(agencies)}) for file key vlbi_obs_vgosdb. Using {_file_vars['agency']}"
                )

        if not "obs_version" in _file_vars and not "acengy" in _file_vars:
            log.fatal(f"No VGOSDB wrapper file found ({config.files.path('vlbi_obs_vgosdb')}).")

    # Sinex file vars
    if "sinex" in config.tech.section_names:
        _file_vars["solution"] = config.tech.sinex.solution.str
        _file_vars["file_agency"] = config.tech.sinex.file_agency.str.lower()

    return _file_vars


@plugins.register_named("log_prefix")
def log_prefix(rundate, session="", **kwargs):
    return f"{pipeline.upper()} {session} {rundate:%Y-%m-%d}"


def _matplotlib_map(dset):
    # Local imports since only this function uses plotting
    import cartopy.crs as ccrs
    from matplotlib import cm
    import matplotlib.pyplot as plt

    stations = dset.unique("station")
    baselines = dset.unique("baseline")
    b_width = {b: dset.num(baseline=b) / dset.num_obs * 40 for b in baselines}
    s_size = [int(dset.num(station=s) / dset.num_obs * 400) for s in stations]

    lat = np.degrees([dset.meta[station]["latitude"] for station in stations])
    lon = np.degrees([dset.meta[station]["longitude"] for station in stations])

    rms = [dset.rms("residual", station=station) for station in stations] if "residual" in dset._fields else None
    plt.figure()
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.set_global()
    ax.coastlines()
    cmap = cm.get_cmap(config.there.colormap.str)

    for baseline in baselines:
        sta_1, _, sta_2 = baseline.partition("/")
        sta_1_lon = np.degrees(dset.meta[sta_1]["longitude"])
        sta_2_lon = np.degrees(dset.meta[sta_2]["longitude"])
        sta_1_lat = np.degrees(dset.meta[sta_1]["latitude"])
        sta_2_lat = np.degrees(dset.meta[sta_2]["latitude"])
        b_rms = dset.rms("residual", baseline=baseline) if "residual" in dset._fields else 0
        color = cmap(b_rms)
        plt.plot(
            [sta_1_lon, sta_2_lon],
            [sta_1_lat, sta_2_lat],
            c=color,
            linestyle="--",
            linewidth=b_width[baseline],
            transform=ccrs.Geodetic(),
            zorder=0,
        )
    plt.scatter(lon, lat, c=rms, transform=ccrs.Geodetic(), cmap=cmap, s=s_size, zorder=10)
    plt.title(f"{dset.meta['input']['session_code']}")
    plt.figtext(0.99, 0.01, "Size: number of observations", horizontalalignment="right")
    if rms is not None:
        cbar = plt.colorbar()
        cbar.set_label("RMS of residual [m]")

    plt.show()


@plugins.register_named("make_map")
def make_map(dset):
    """Plots station positions on a global map."""

    web_map = config.where.get(
        "web_map",
        section=pipeline,
        value=util.read_option_value("--web_map", default=None),  # TODO: add this to mg_config
        default=False,
    ).bool

    if web_map:
        import webbrowser  # Local import (Only needed if this function is called)

        map_path = config.files.path("output_web_map", file_vars={**dset.vars, **dset.analysis})
        if not map_path.exists():
            writers.write_one("web_map", dset=dset)
        webbrowser.open(map_path.as_uri())

    plot_map = config.where.get(
        "plot_map",
        section=pipeline,
        value=util.read_option_value("--plot_map", default=None),  # TODO: add this to mg_config
        default=False,
    ).bool

    if plot_map:
        _matplotlib_map(dset)


### Stages


@plugins.register
def read(stage, dset):
    """Read VLBI data

    Args:
        dset (Dataset):  Dataset
    """
    obs.get(dset)
    log.info(f"Parsed {dset.num_obs} observations")
    dset.write_as(stage=stage, label=0)


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

    dset.write_as(stage=stage, label=0)


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
    site.calculate_site("site", dset)
    delta_pos = site.add("site", dset)

    dset.site_pos_1[:] = (dset.site_pos_1.gcrs + delta_pos[0].gcrs).trs
    dset.site_pos_2[:] = (dset.site_pos_2.gcrs + delta_pos[1].gcrs).trs
    log.blank()

    # Run models for each term of the observation equation
    log.info(f"Calculating theoretical delays")
    delay.calculate_delay("delay", dset)

    delta_delay = delay.add("delay", dset)
    dset.add_float("obs", val=dset.observed_delay, unit="meter", write_level="operational")
    dset.add_float("calc", val=delta_delay, unit="meter", write_level="operational")
    dset.add_float("residual", val=dset.obs - dset.calc, unit="meter", write_level="operational")
    log.blank()

    # Estimate clock polynomial
    log.info(f"Calculating clock polynomials")
    max_iterations = config.tech.calculate_max_iterations.int
    outlier_limit = config.tech.calculate_outlier_limit.float
    store_outliers = config.tech.store_outliers.bool

    for iter_num in itertools.count(start=1):
        if dset.num_obs == 0:
            break
        delay.calculate_delay("delay_corr", dset)
        delta_correction = delay.add("delay_corr", dset)

        dset.calc[:] = dset.calc + delta_correction
        dset.residual[:] = dset.obs - dset.calc
        rms = dset.rms("residual")
        log.info(f"{dset.num_obs} observations, residual = {rms:.4f}")

        # Store results
        dset.write_as(stage=stage, label=iter_num - 1)

        # Detect and remove extreme outliers
        idx = np.abs(dset.residual) < outlier_limit * rms
        if iter_num > max_iterations or idx.all():
            break

        if store_outliers:
            bad_idx = np.logical_not(idx)
            log.info(f"Adding {np.sum(bad_idx)} observations to ignore_observation")
            bad_obs = np.char.add(np.char.add(dset.time.utc.iso[bad_idx], " "), dset.baseline[bad_idx]).tolist()
            with config.update_tech_config(dset.analysis["rundate"], pipeline, session=dset.vars["session"]) as cfg:
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
    max_iterations = config.tech.estimate_max_iterations.int

    for iter_num in itertools.count(start=1):
        partial_vectors = estimation.partial_vectors(dset, "estimate_method")
        obs_noise = dset.observed_delay_ferr ** 2 + np.nan_to_num(dset.iono_delay_ferr) ** 2 + 0.01 ** 2
        log.info(
            f"Estimating parameters for iteration {iter_num} using Kalman Filter and continuous piecewise linear functions"
        )
        estimation.call("estimate_method", dset=dset, partial_vectors=partial_vectors, obs_noise=obs_noise)
        rms = dset.rms("residual")
        log.info(f"{dset.num_obs} observations, postfit residual = {rms:.4f}")
        dset.write_as(stage=stage, label=iter_num - 1)
        if iter_num >= max_iterations:
            break

        # Detect and remove outliers
        num_obs_before = dset.num_obs
        independent = config.tech.estimate_obs_rejectors_independent.bool
        dset = estimation.apply_observation_rejectors("estimate_obs_rejectors", dset, independent)
        log.blank()
        if dset.num_obs == num_obs_before or dset.num_obs == 0:
            break

    log.blank()
    estimation.solve_neq(dset)
    dset.write()


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
