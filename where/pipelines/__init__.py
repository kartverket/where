"""Frameworks for the available pipelines

Description:
------------

To add a new pipeline, simply create a new .py-file in the pipelines-directory. Each pipeline is split into several
parts called stages. These stages will be called in the order they are registered in the source file (by using the
:func:`register` decorator).

The idea is that each stage should be a self-contained step in the analysis. Typically data are transferred between
stages by being stored to and read from disk. Although this adds a little overhead, it has the advantage that those
data are available after the analysis is finished.

"""
# Standard library import
import sys

# Where imports
from where import reports
from where import setup
from where.lib import cache
from where.lib import config
from where.lib import dependencies
from where.lib import exceptions
from where.lib import files
from where.lib import log
from where.lib import plugins
from where.lib.timer import timer
from where.lib import util


@cache.function
def options():
    """List the command line options for starting the different pipelines

    Returns:
        Dict:  Command line options pointing to pipelines
    """
    options = dict()
    plugin_files = plugins.list_all(package_name=__name__)

    for pipeline in plugin_files:
        try:
            pipeline_options = plugins.call_one(package_name=__name__, plugin_name=pipeline, part="options")
        except exceptions.UnknownPluginError:
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
    return plugins.list_parts(package_name=__name__, plugin_name=pipeline)


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
    if pipeline is not None:  # Can not use plugins.doc because we want the module doc string
        return sys.modules[f"{__name__}.{pipeline}"].__doc__.split("\n")[0]

    doc_list = list()
    all_options = options()

    for pipeline in names():
        pipeline_opts = ", ".join(o for o, p in all_options.items() if p == pipeline)
        doc_list.append(f"{pipeline_opts:21s}{{action}} {doc(pipeline)}")

    # Add a custom formatter to customize the action
    class DocString(str):

        def __format__(self, format_spec):
            if not format_spec:
                format_spec = "Run a"
            return self.format(action=format_spec)

    return DocString("\n".join(doc_list))


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
        return plugins.call_one(
            package_name=__name__, plugin_name=pipeline, part="validate_session", rundate=rundate, session=session
        )
    except exceptions.UnknownPluginError:
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
        return plugins.call_one(package_name=__name__, plugin_name=pipeline, part="list_sessions", rundate=rundate)
    except exceptions.UnknownPluginError:
        return [""]  # If sessions is not defined in the pipeline, return a list with one unnamed session


def file_vars():
    """Get a list of file variables for the current pipeline

    The active analysis variables are also made available, but may be overridden by the pipeline.
    """
    file_vars = dict(config.analysis.config.as_dict(), **config.date_vars(config.analysis.rundate.date))
    pipeline_file_vars = plugins.call_one(
        package_name=__name__, plugin_name=config.analysis.tech.str, part="file_vars"
    )
    file_vars.update(pipeline_file_vars)

    return file_vars


def call(pipeline_, stage_, **stage_args):
    """Call one stage of the pipeline

    The underscore-postfix of `pipeline_` and `stage_` is used so that it does not interfere with the stage_args.

    Args:
        pipeline_ (String):  The pipeline.
        stage_ (String):     The stage.
        stage_args:          Arguments that will be passed to the stage-function.

    Returns:
        The return value of the stage-function.
    """
    return plugins.call_one(package_name=__name__, plugin_name=pipeline_, part=stage_, **stage_args)


def run(rundate, pipeline, session=""):
    """Run a Where pipeline for a given date and session

    Args:
        rundate:   Rundate of analysis.
        pipeline:  Pipeline used for analysis.
        session:   Session in analysis.
    """
    if not setup.has_config(rundate, pipeline, session):
        log.fatal(f"No configuration found for {pipeline.upper()} {session} {rundate.strftime(config.FMT_date)}")

    # Set up tech config and file logging
    config.init(rundate=rundate, tech_name=pipeline, session=session)
    log.file_init(log_path=files.path("log"))

    # Read which stages to skip from technique configuration file.
    skip_stages = config.tech.get("skip_stages", default="").list

    # Register filekey suffix
    filekey_suffix = config.tech.filekey_suffix.list
    if filekey_suffix:
        files.use_filelist_profiles(*filekey_suffix)

    # Find which stages we will run analysis for
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
    dep_fast = config.where.files.dependencies_fast.bool
    for prev_stage, stage in zip([None] + stage_list, stage_list):

        # Skip stages where no dependencies have changed
        if not (
            dependencies.changed(fast_check=dep_fast, rundate=rundate, tech=pipeline, session=session, stage=stage)
            or util.check_options("-F", "--force")
        ):
            log.info(f"Not necessary to run {stage} for {pipeline.upper()} {rundate.strftime(config.FMT_date)}")
            continue

        # Report on the stage
        reports.report.start_section(stage)
        reports.report.text("header", stage.replace("_", " ").title())
        if prev_stage:
            log.blank()  # Empty line for visual clarity

        # Set up dependencies. Add dependencies to previous stage and config file
        dependencies.init(fast_check=dep_fast, session=session, stage=stage)
        dependencies.add(files.path("model_run_depends", file_vars=dict(session=session, stage=prev_stage)))
        dependencies.add(*config.tech.sources)

        # Call the current stage. Skip rest of stages if current stage returns False (compare with is since by
        # default stages return None)
        do_next_stage = call(
            pipeline, stage, rundate=rundate, session=session, prev_stage=prev_stage, stage=stage, logger=log.info
        )
        dependencies.write()
        if do_next_stage is False:
            break  # TODO, this does not work together with dependencies changed ...

    # Publish files for session
    files.publish_files()
    session_timer.end()

    # Store configuration to library
    setup.store_config_to_library(rundate, pipeline, session)

    # Write reports specified in config
    reports.write(rundate, pipeline)

    # Write requirements to file for reproducibility
    util.write_requirements()
