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

# Standard library imports
from functools import lru_cache

# Midgard imports
from midgard.dev import plugins
from midgard.dev import exceptions as mg_exceptions
from midgard.files import dependencies

# Where imports
from where import setup
from where.data import dataset3 as dataset
from where.lib import config
from where.lib import exceptions
from where.lib import log
from where.lib import util


@lru_cache()
def options():
    """List the command line options for starting the different pipelines

    Returns:
        Dict:  Command line options pointing to pipelines
    """
    opts = dict()
    plugin_files = plugins.names(package_name=__name__)

    for pipeline in plugin_files:
        try:
            pipeline_options = plugins.call(package_name=__name__, plugin_name=pipeline, part="options")
        except mg_exceptions.UnknownPluginError:
            continue
        opts.update({opt: pipeline for opt in pipeline_options})

    return opts


@lru_cache()
def names():
    """List the names of the pipelines available

    Returns:
        List: Strings with the names of the available pipelines.
    """
    return sorted({pipeline for pipeline in options().values()})


@lru_cache()
def stages(pipeline):
    """Get a list of stages for a given pipeline

    Args:
        pipeline (String): The pipeline.

    Returns:
        Iterator: Strings with the names of the stages.
    """
    return plugins.parts(package_name=__name__, plugin_name=pipeline)


@lru_cache()
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


def get_args(rundate, pipeline, input_args=None):
    """Get a list of sessions for a given rundate for a pipeline

    Args:
        rundate (Date):     The model run date.
        pipeline (String):  Name of pipeline.
        input_args:         List of arguments from the command line

    Returns:
        List: List of command line arguments
    """
    try:
        return plugins.call(
            package_name=__name__, plugin_name=pipeline, part="get_args", rundate=rundate, input_args=input_args
        )
    except mg_exceptions.UnknownPluginError:
        log.warn(f"Pipeline {pipeline} has not defined function get_args")
        return input


def file_vars():
    """Get a list of file variables for the current pipeline

    The active analysis variables are also made available, but may be overridden by the pipeline.
    """
    _file_vars = dict(config.analysis.config.as_dict(), **config.date_vars(config.analysis.rundate.date))
    pipeline_file_vars = plugins.call(
        package_name=__name__, plugin_name=config.analysis.pipeline.str, part="file_vars", file_vars=_file_vars
    )
    _file_vars.update(pipeline_file_vars)

    return _file_vars


def make_map(dset):
    """Make and show a basic matplotlib plot relevant for the current pipeline"""
    try:
        plugins.call(package_name=__name__, plugin_name=dset.vars["pipeline"], part="make_map", dset=dset)
    except mg_exceptions.UnknownPluginError:
        log.warn(f"Pipeline {dset.vars['pipeline']} has not defined function make_map")


def paths(label_pattern, pipeline=None):
    """Get a list of dependent file paths with a given label

    Args:
        label_pattern:  String with label or regular expression (e.g. 'gnss_rinex_nav_[MGE]' or 'gnss_rinex_nav_.').
        pipeline:       Pipeline used for analysis.

    Returns:
        Set:  Set of file paths.
    """
    pipeline = config.analysis.pipeline.str if pipeline is None else pipeline
    file_paths = list()
    for stage in stages(pipeline):
        dep_path = config.files.path("depends", file_vars=dict(stage=stage))
        file_paths.extend(dependencies.get_paths_with_label(dep_path, label_pattern))

    return set(file_paths)


def run(rundate, pipeline, *args, **kwargs):
    """Run a Where pipeline for a given date and session

    Args:
        rundate:   Rundate of analysis.
        pipeline:  Pipeline used for analysis.
        session:   Session in analysis.
    """

    if not setup.has_config(rundate, pipeline, *args, **kwargs):
        log.fatal(f"No configuration found for {pipeline.upper()} {rundate.strftime(config.FMT_date)}")

    # Set up config
    config.init(rundate, pipeline, **kwargs)

    # Register filekey suffix
    filekey_suffix = config.tech.filekey_suffix.list
    if filekey_suffix:
        config.files.profiles = filekey_suffix

    # Validate input arguments
    try:
        prefix = plugins.call(
            package_name=__name__, plugin_name=pipeline, part="validate_args", rundate=rundate, **kwargs
        )
    except mg_exceptions.UnknownPluginError:
        log.warn(f"Pipeline {pipeline} has not defined function 'validate_args'")
    except exceptions.InvalidArgsError as err:

        from where.tools import delete

        # Clean up {placeholder} directories created by config
        delete.delete_analysis(rundate, pipeline, **kwargs)
        log.fatal(err)

    # Set up console logger and start file logger
    try:
        prefix = plugins.call(
            package_name=__name__, plugin_name=pipeline, part="log_prefix", rundate=rundate, **kwargs
        )
    except mg_exceptions.UnknownPluginError:
        log.warn(f"Pipeline {pipeline} has not defined function 'log_prefix'")
        prefix = ""

    log_cfg = config.where.log
    log.init(log_level=log_cfg.default_level.str, prefix=prefix)
    if log_cfg.log_to_file.bool:
        log.file_init(
            file_path=config.files.path("log"),
            log_level=log_cfg.default_level.str,
            prefix=prefix,
            rotation=log_cfg.number_of_log_backups.int,
        )

    # Update analysis config and file variables
    config.set_analysis(rundate, pipeline=pipeline, **kwargs)
    config.set_file_vars(file_vars())

    log.blank()  # Empty line for visual clarity

    # Read which stages that should be executed once for each iterable
    skip_stages = config.tech.skip_stages.list
    stage_iterate = config.tech.stage_iterate.list
    dset_list = []
    dset = None

    if stage_iterate:
        # Read which list should be iterated over and the placeholder name of each entry
        iterate_over, _, var_name = config.tech.stage_iterate_over.str.partition(":")
        var_name = var_name.strip()

        # Iterate
        for item in config.tech[iterate_over].list:
            kwargs[var_name] = item
            log.blank()
            log.info(f"***** Running {item} *****")

            for prev_stage, stage in zip([None] + stage_iterate, stage_iterate):
                if stage not in skip_stages:
                    dset = run_stage(rundate, pipeline, dset, stage, prev_stage, **kwargs)

            if dset is not None:
                dset_list.append(dset)
                dset = None
        kwargs[var_name] = "combined"

    if dset_list:
        dset_list[0].merge_with(*dset_list[1:], sort_by="time")
        dset = dset_list[0]
        if len(dset_list) > 1:
            log.info(f"Combining dataset for {len(dset_list)} {iterate_over}")
            dset.write_as(stage=stage_iterate[-1], label=2, **kwargs)

    # Read which stages that should be executed once
    stage_once = config.tech.stage_once.list
    # Find which stages we will run analysis for
    if not stage_once and not stage_iterate:
        stage_list = [s for s in stages(pipeline)]
        prev_stage_start = None
    else:
        stage_list = [s for s in stage_once]
        prev_stage_start = stage_iterate[-1] if stage_iterate else None

    for prev_stage, stage in zip([prev_stage_start] + stage_list, stage_list):
        if stage not in skip_stages:
            dset = run_stage(rundate, pipeline, dset, stage, prev_stage, **kwargs)
            log.blank()

        if dset is not None and dset.num_obs == 0:
            log.warn(f"No observations in dataset after {stage} stage.")
            break

    # Store configuration to library
    setup.store_config_to_library(rundate, pipeline, **kwargs)

    # Write requirements to file for reproducibility
    util.write_requirements()


def run_stage(rundate, pipeline, dset, stage, prev_stage, **kwargs):
    # Skip stages where no dependencies have changed
    dep_path = config.files.path("depends", file_vars={**kwargs, "stage": stage})
    if not (dependencies.changed(dep_path) or util.check_options("-F", "--force")):
        log.info(f"Not necessary to run {stage} for {pipeline.upper()} {rundate.strftime(config.FMT_date)}")
        return

    if dset is None:
        try:
            # Read dataset from disk if it exists
            dset = dataset.Dataset.read(rundate=rundate, pipeline=pipeline, stage=prev_stage, label="last", **kwargs)
        except (OSError, ValueError):
            # Create emtpy dataset
            dset = dataset.Dataset(rundate=rundate, pipeline=pipeline, **kwargs)

    # Set up dependencies. Add dependencies to previous stage and config file
    dependencies.init(dep_path)
    if prev_stage is not None:
        dependencies.add(config.files.path("depends", file_vars={**kwargs, "stage": prev_stage}), label="depends")
    dependencies.add(*config.tech.sources, label="config")
    # Delete old datasets for this stage
    dset.delete_stage(stage, **kwargs)

    # Call the current stage. Skip rest of stages if current stage returns False (compare with is since by
    # default stages return None)
    plugins.call(
        package_name=__name__, plugin_name=pipeline, part=stage, stage=stage, dset=dset, plugin_logger=log.info
    )
    dependencies.write()

    return dset
