"""Frameworks for the available pipelines

Description:
------------

To add a new pipeline, simply create a new .py-file in the pipelines-directory. Each pipeline is split into several
parts called stages. These stages will be called in the order they are registered in the source file (by using the
`@plugins.register` decorator).

The idea is that each stage should be a self-contained step in the analysis. Typically data are transferred between
stages by being stored to and read from disk. Although this adds a little overhead, it has the advantage that those
data are available after the analysis is finished.
"""

# Midgard imports
from midgard.dev import plugins
from midgard.dev import exceptions as mg_exceptions
from midgard.files import dependencies

# Where imports
from where import reports
from where import setup
from where.data import dataset
from where.lib import cache
from where.lib import config
from where.lib import exceptions
from where.lib import files
from where.lib import log
from where.lib.timer import timer
from where.lib import util


@cache.function
def options():
    """List the command line options for starting the different pipelines

    Returns:
        Dict:  Command line options pointing to pipelines
    """
    options = dict()
    plugin_files = plugins.names(package_name=__name__)

    for pipeline in plugin_files:
        try:
            pipeline_options = plugins.call(package_name=__name__, plugin_name=pipeline, part="options")
        except mg_exceptions.UnknownPluginError:
            continue
        options.update({opt: pipeline for opt in pipeline_options})

    return options


@cache.function
def names():
    """List the names of the pipelines available

    Returns:
        List: Strings with the names of the available pipelines.
    """
    return sorted({pipeline for pipeline in options().values()})


@cache.function
def stages(pipeline):
    """Get a list of stages for a given pipeline

    Args:
        pipeline (String): The pipeline.

    Returns:
        Iterator: Strings with the names of the stages.
    """
    return plugins.parts(package_name=__name__, plugin_name=pipeline)


@cache.function
def doc(pipeline=None):
    """Get doc for one or all pipelines

    If pipeline is specified, return the doc-string for that pipeline, otherwise return a string with short
    documentation for all available pipelines.

    Args:
        pipeline (String):  The pipeline to document (optional).

    Returns:
        String:  Documentation of one or all pipelines.
    """
    # Add a custom formatter to customize the action
    class DocString(str):
        def __format__(self, format_spec):
            if not format_spec:
                format_spec = "Run"
            return self.format(action=format_spec)

    # Single pipeline
    if pipeline is not None:
        doc_str = plugins.doc(package_name=__name__, plugin_name=pipeline, long_doc=False, use_module=True)
        return DocString(f"{{action}} {doc_str}")

    # All pipelines
    doc_list = list()
    all_options = options()

    for pipeline in names():
        pipeline_opts = ", ".join(o for o, p in all_options.items() if p == pipeline)
        doc_list.append(f"{pipeline_opts:<20s} {doc(pipeline):__action__}")

    return DocString("\n".join(doc_list).replace("__action__", "{action}"))


def get_from_options():
    """Read pipeline from command line options

    Uses the `options`-plugin in each pipeline to find which pipeline is called. Raises UnknownPipelineError if there
    are no recognized pipelines in the command line options.

    Returns:
        String:  Name of pipeline.
    """
    available_pipelines = options()
    pipeline_option = util.check_options(*available_pipelines)
    if not pipeline_option:
        raise exceptions.UnknownPipelineError(
            f"No pipeline specified.\nUse one of {', '.join(available_pipelines)} to specify the pipeline"
        )

    return available_pipelines[pipeline_option]


def get_session(rundate, pipeline):
    """Read session from command line options

    The session is validated for the given pipeline. Uses the `validate_session`-plugin for validation.

    Args:
        pipeline (String):  Name of pipeline.

    Returns:
        String:  Name of session.
    """
    session = util.read_option_value("--session", default="")
    try:
        return plugins.call(
            package_name=__name__, plugin_name=pipeline, part="validate_session", rundate=rundate, session=session
        )
    except mg_exceptions.UnknownPluginError:
        return session  # Simply return session if it can not be validated


@cache.function
def list_sessions(rundate, pipeline):
    """Get a list of sessions for a given rundate for a pipeline

    Args:
        rundate (Date):     The model run date.
        pipeline (String):  Name of pipeline.

    Returns:
        List: Strings with the names of the sessions.
    """
    try:
        return plugins.call(package_name=__name__, plugin_name=pipeline, part="list_sessions", rundate=rundate)
    except mg_exceptions.UnknownPluginError:
        return [""]  # If sessions is not defined in the pipeline, return a list with one unnamed session


def file_vars():
    """Get a list of file variables for the current pipeline

    The active analysis variables are also made available, but may be overridden by the pipeline.
    """
    file_vars = dict(config.analysis.config.as_dict(), **config.date_vars(config.analysis.rundate.date))
    pipeline_file_vars = plugins.call(package_name=__name__, plugin_name=config.analysis.tech.str, part="file_vars")
    file_vars.update(pipeline_file_vars)

    return file_vars


def paths(label_pattern, pipeline=None):
    """Get a list of dependent file paths with a given label

    Args:
        label_pattern:  String with label or regular expression (e.g. 'gnss_rinex_nav_[MGE]' or 'gnss_rinex_nav_.').
        pipeline:       Pipeline used for analysis.

    Returns:
        Set:  Set of file paths.
    """
    pipeline = config.analysis.tech.str if pipeline is None else pipeline
    paths = list()
    for stage in stages(pipeline):
        dep_path = files.path("depends", file_vars=dict(stage=stage))
        paths.extend(dependencies.get_paths_with_label(dep_path, label_pattern))

    return set(paths)


def run(rundate, pipeline, session=""):
    """Run a Where pipeline for a given date and session

    Args:
        rundate:   Rundate of analysis.
        pipeline:  Pipeline used for analysis.
        session:   Session in analysis.
    """
    if not setup.has_config(rundate, pipeline, session):
        log.fatal(f"No configuration found for {pipeline.upper()} {session} {rundate.strftime(config.FMT_date)}")

    # Set up session config
    config.init(rundate=rundate, tech_name=pipeline, session=session)

    # Set up prefix for console logger and start file logger
    log_cfg = config.where.log
    prefix = f"{pipeline.upper()} {session} {rundate:%Y-%m-%d}"
    log.init(log_level=log_cfg.default_level.str, prefix=prefix)
    if log_cfg.log_to_file.bool:
        log.file_init(
            file_path=files.path("log"),
            log_level=log_cfg.default_level.str,
            prefix=prefix,
            rotation=log_cfg.number_of_log_backups.int,
        )

    # Read which stages to skip from technique configuration file.
    skip_stages = config.tech.get("skip_stages", default="").list

    # Register filekey suffix
    filekey_suffix = config.tech.filekey_suffix.list
    if filekey_suffix:
        config.files.profiles = filekey_suffix

    # Find which stages we will run analysis for
    # TODO: Specify stage_list in config
    stage_list = [s for s in stages(pipeline) if s not in skip_stages]

    # Start file logging and reporting
    reports.report.init(sessions=[session])
    reports.report.start_session(session)
    reports.report.text("header", session.replace("_", " ").title())

    # Update analysis config and file variables
    config.set_analysis(rundate=rundate, tech=pipeline, analysis=pipeline, session=session)
    config.set_file_vars(file_vars())

    # Log the name of the session
    log.blank()  # Empty line for visual clarity
    log.info(f"Start session {session}")
    session_timer = timer(f"Finish session {session} in")
    session_timer.start()

    # Run stages, keep track of previous stage
    dset = None
    dep_fast = config.where.files.dependencies_fast.bool
    for prev_stage, stage in zip([None] + stage_list, stage_list):

        # Skip stages where no dependencies have changed
        dep_path = files.path("depends", file_vars=dict(stage=stage))
        if not (dependencies.changed(dep_path, fast_check=dep_fast) or util.check_options("-F", "--force")):
            log.info(f"Not necessary to run {stage} for {pipeline.upper()} {rundate.strftime(config.FMT_date)}")
            continue
        elif dset is None:
            # Create or read dataset
            empty = stage == stage_list[0]
            dset = dataset.Dataset(
                rundate, tech=pipeline, stage=prev_stage, dataset_name=session, dataset_id="last", empty=empty
            )

        # Report on the stage
        reports.report.start_section(stage)
        reports.report.text("header", stage.replace("_", " ").title())
        if prev_stage:
            log.blank()  # Empty line for visual clarity

        # Set up dependencies. Add dependencies to previous stage and config file
        dependencies.init(dep_path, fast_check=dep_fast)
        dependencies.add(files.path("depends", file_vars=dict(stage=prev_stage)), label="depends")
        dependencies.add(*config.tech.sources, label="config")

        # Delete old datasets for this stage
        dset.delete_from_file(stage=stage, dataset_id="all")

        # Call the current stage. Skip rest of stages if current stage returns False (compare with is since by
        # default stages return None)
        plugins.call(
            package_name=__name__, plugin_name=pipeline, part=stage, stage=stage, dset=dset, plugin_logger=log.info
        )
        dependencies.write()
        if dset.num_obs == 0:
            log.warn(f"No observations in dataset after {stage} stage. Exiting pipeline")
            break
    else:  # Only done if loop does not break (all stages finish normally)
        # Publish files for session
        files.publish_files()

    session_timer.end()

    # Store configuration to library
    setup.store_config_to_library(rundate, pipeline, session)

    # Write reports specified in config
    reports.write(rundate, pipeline)

    # Write requirements to file for reproducibility
    util.write_requirements()
