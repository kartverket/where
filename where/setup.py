#!/usr/bin/env python3
"""Set up the Where program to do analysis of space geodetic data

Usage::

    {exe:setup} date pipeline [--session=session] [options]

The program requires a date. Typically, the date is given in the format
`<year month day>` (for example 2015 8 4). However, it is also possible to
specify the date as `<year day-of-year>` (for example 2015 216) by also adding
the option `--doy`.
In addition, one pipeline must be specified. See below for available pipelines.

===================  ===========================================================
Pipeline             Description
===================  ===========================================================
{pipelines_doc:Set up}
===================  ===========================================================

Furthermore, the following options are recognized:

===================  ===========================================================
Option               Description
===================  ===========================================================
-D, --delete         Delete existing analysis results.
    --doy            Specify date as <year day-of-year>.
-E, --edit           Edit the configuration of an analysis.
-N, --new            Start a new analysis (combine with -A or -D).
-T, --showtb         Show traceback if the program crashes.
--id=analysisid      Add a special analysis id (to run several versions of the
                     same analysis simultaneously).
--profile=name       Use config settings specfied for a given profile, for
                     instance --profile=vascc for a VLBI analysis.
--session=session    Set up analysis for the given session.
--user=username      Run as username. Does not need to be an existing username
                     on the system.
--debug, ...         Show additional debug information. Other flags such as
                     --all, --debug, --time, --dev, --info, --out, --warn,
                     --check, --error, --fatal, --none are also allowed, and
                     shows differing amounts of information as the program runs.
--version            Show version information and exit.
-h, --help           Show this help message and exit.
===================  ===========================================================


Description:
------------

This program is used to set up a Where analysis. The program will create a
configuration file for the model run and display it. This configuration can be
changed, either by editing the given configuration file or by rerunning
{exe:setup} with different configuration options.

See :doc:`user_guide_where` for more information.


Examples:
---------

Here are some concrete examples of how to run the program:

Set up a VLBI analysis for August 4 2015::

    {exe:setup} 2015 8 4 -v

Set up an SLR analysis for September 1 2015 using day-of-year::

    {exe:setup} 2015 242 --slr --doy

Change the spam option of the GNSS analysis::

    {exe:setup} 2016 3 1 -g --spam=ham


Current Maintainers:
--------------------

{maintainers}

Version: {version}

"""
# Standard library imports
from datetime import datetime
import sys

# Midgard imports
from midgard.config import config as mg_config
from midgard.dev import console

# Where imports
from where.tools import delete
from where import pipelines
from where.lib import config
from where.lib import files
from where.lib import log
from where.lib.timer import timer
from where.lib import util

# Optional import of editor and seaborn, add dummy methods in case they are not installed
from midgard.dev import optional

editor = optional.optional_import("editor")


@timer(f"Finish {util.get_program_name()} in")
@util.no_traceback
def main():
    """Parse command line options and set up an Where analysis

    Do simple parsing of command line arguments. Set up config-files and show the configuration.
    """
    # Start logging
    log.init()

    # Read command line options
    if util.check_options("--doy"):
        rundate = util.parse_args("doy", doc_module=__name__)
    else:
        rundate = util.parse_args("date", doc_module=__name__)
    pipeline = pipelines.get_from_options()
    session = pipelines.get_session(rundate, pipeline)

    # Set up the configuration for the analysis
    setup_config(rundate, pipeline, session)

    # Show current configuration
    show_config(rundate, pipeline, session)

    # Store configuration in library
    store_config_to_library(rundate, pipeline, session)


def setup_config(rundate, pipeline, session):
    """Set up configuration for a Where analysis

    """
    # Set the correct profile
    profile = util.read_option_value("--profile", default="")
    config.where.profiles = profile.split() + [pipeline]

    # Should a new analysis be started?
    start_new = util.check_options("-N", "--new")

    # Delete an analysis
    if util.check_options("-D", "--delete"):
        delete.delete_analysis(rundate, pipeline, session)
        if not start_new:
            raise SystemExit

    # Create configuration of a new analysis
    if start_new or not has_config(rundate, pipeline, session):
        create_config(rundate, pipeline, session)
    elif util.check_options("--profile"):  # Warning if --profile option is ignored
        profile_opt = f"--profile={util.read_option_value('--profile', default='')}"
        log.warn(f"Configuration already exists, option '{profile_opt}' ignored")

    # Update configuration based on command line options
    update_config(rundate, pipeline, session)

    # Edit configuration manually
    if util.check_options("-E", "--edit"):
        edit_config(rundate, pipeline, session)

    # Show current configuration
    if util.check_options("-S", "--show-config"):
        show_config(rundate, pipeline, session)
        raise SystemExit


def has_config(rundate, pipeline, session):
    """Test whether the configuration of a Where analysis exists

    """
    return _config_path(rundate, pipeline, session).exists()


def create_config(rundate, pipeline, session):
    """Create the configuration of a Where analysis

    """
    # Create a new configuration and copy all and pipeline sections
    cfg_path = _config_path(rundate, pipeline, session)
    cfg = mg_config.Configuration(pipeline)
    cfg.update_from_config_section(config.where.all, section=pipeline)
    cfg.update_from_config_section(config.where[pipeline], section=pipeline)
    log.info(f"Creating new configuration at '{cfg_path}' based on {', '.join(cfg.sources)}")
    cfg.write_to_file(cfg_path, metadata=False)

    # Update configuration settings from library
    for section in read_from_library(rundate, pipeline, session):
        cfg.update_from_config_section(section, section.name)

    # Write updated configuration to file
    cfg.write_to_file(cfg_path, metadata=False)

    # Add timestamp and creation note
    add_timestamp(rundate, pipeline, session, "created")

    # Add new dependent sections from newly created config
    add_sections(rundate, pipeline, session)


def update_config(rundate, pipeline, session):
    """Update the configuration of a Where analysis

    """
    cfg_path = _config_path(rundate, pipeline, session)
    ts_before = files.get_timestamp(cfg_path)

    # Update with command line options
    with mg_config.Configuration.update_on_file(_config_path(rundate, pipeline, session)) as cfg:
        cfg.master_section = pipeline
        cfg.update_from_options(_clean_sys_argv())

    # Add timestamp and updated note
    if files.get_timestamp(cfg_path) != ts_before:
        add_timestamp(rundate, pipeline, session, "last update")

    # Add new dependent sections from command line options
    add_sections(rundate, pipeline, session)


def edit_config(rundate, pipeline, session):
    """Update the configuration of a Where analysis

    """
    cfg_path = _config_path(rundate, pipeline, session)
    ts_before = files.get_timestamp(cfg_path)

    # Open config file in an editor
    editor.edit(str(cfg_path))

    if files.get_timestamp(cfg_path) != ts_before:
        # Add timestamp and edited note
        add_timestamp(rundate, pipeline, session, "last update")

        # Add new dependent sections from manual edit
        add_sections(rundate, pipeline, session)


def add_sections(rundate, pipeline, session):
    """Update the configuration with sections with settings for models, cleaners etc

    Todo: Figure out how to work with metadata across profiles

    """
    cfg_path = _config_path(rundate, pipeline, session)
    ts_before = files.get_timestamp(cfg_path)

    # Add dependent sections that are not already included
    with mg_config.Configuration.update_on_file(cfg_path, metadata=False) as cfg:
        for section in _dependent_sections(cfg[pipeline]):
            if section.name not in cfg.section_names:
                cfg.update_from_config_section(section)
            else:
                for key, entry in section.items():
                    if key not in cfg[section.name]:
                        cfg.update(section.name, key, entry.str, source=entry.source, meta=entry.meta)

    # Add timestamp and updated note
    if files.get_timestamp(cfg_path) != ts_before:
        add_timestamp(rundate, pipeline, session, "last update")


def show_config(rundate, pipeline, session):
    """Show the configuration of a Where analysis

    """
    line = "=" * console.columns()

    # Warn about missing session
    if not has_config(rundate, pipeline, session):
        log.warn(f"No configuration found for {pipeline.upper()} {session} {rundate.strftime(config.FMT_date)}")

    # Read configuration from file
    else:
        cfg = _read_config(rundate, pipeline, session)

        # Print configuration to console
        print(line)
        print(f"{pipeline.upper()} {session} {rundate.strftime(config.FMT_date)}\n")
        print(cfg)
        print(f"\nConfig file at {', '.join(cfg.sources)}")

    # Add instructions about how to update configuration
    print(line)
    pipeline_opt = [o for o, p in pipelines.options().items() if p == pipeline and o.startswith("--")][0]
    cmd = f"{util.get_program_name()} {rundate.year} {rundate.month} {rundate.day} {pipeline_opt} --session={session}"
    print(f"Use '{cmd} --edit' to edit configuration manually")
    print(f"    '{cmd} --<key>=<value>' to update an entry in the [{pipeline}] section")
    print(f"    '{cmd} --<section>:<key>=<value>' to update an entry in a specific section")


def add_timestamp(rundate, pipeline, session, timestamp_key):
    """Write or update a timestamp to file

    Args:
        rundate:        Rundate of analysis.
        pipeline:       Pipeline used for analysis.
        session:        Session in analysis.
        timestamp_key:  Key denoting timestamp.
    """
    # Find timestamp file
    file_vars = config.create_file_vars(rundate, pipeline, session)
    ts_path = files.path("timestamp", file_vars=file_vars)

    # Add timestamp with update note to timestamp file
    with mg_config.Configuration.update_on_file(ts_path) as ts_cfg:
        timestamp = f"{datetime.now().strftime(config.FMT_datetime)} by {util.get_program_info()}"
        ts_cfg.update("timestamps", timestamp_key, timestamp, source=__file__)


def read_from_library(rundate, pipeline, session):
    cfg = _read_config(rundate, pipeline, session)
    cfg.update_from_options(allow_new=True)
    if not cfg.read_from_library.bool:
        raise StopIteration

    file_vars = config.create_file_vars(rundate, pipeline, session)
    lib_path = files.path("config_library", file_vars=file_vars)
    lib_cfg = mg_config.Configuration.read_from_file("library", lib_path)

    for section in lib_cfg.sections:
        yield section


def store_config_to_library(rundate, pipeline, session):
    cfg = _read_config(rundate, pipeline, session)
    if not cfg.write_to_library.bool:
        return

    file_vars = config.create_file_vars(rundate, pipeline, session)
    lib_path = files.path("config_library", file_vars=file_vars)
    lib_cfg = mg_config.Configuration("library")

    for section in cfg.sections:
        for key, entry in section.items():  # Todo: Make ConfigurationSection iterable
            if "library" in entry.meta or "library" in config.where.get(key, section=section.name, default="").meta:
                lib_cfg.update(section.name, key, entry.str, source=entry.source)
                # Todo: Only store entries different from default (issue: profiles?)

    lib_cfg.write_to_file(lib_path)


def _read_config(rundate, pipeline, session):
    """Read the configuration of a Where analysis from file

    Todo: Add this as a classmethod on Configuration

    Args:
        rundate:   Rundate of analysis.
        pipeline:  Pipeline used for analysis.
        session:   Session in analysis.

    Returns:
        Configuration of Where analysis.
    """
    if not has_config(rundate, pipeline, session):
        raise FileNotFoundError(
            f"No configuration found for {pipeline.upper()} {session} {rundate.strftime(config.FMT_date)}"
        )

    cfg = mg_config.Configuration.read_from_file(pipeline, _config_path(rundate, pipeline, session))
    cfg.master_section = pipeline

    return cfg


def _config_path(rundate, pipeline, session):
    """The path to the configuration of a Where analysis

    Todo: Move this to lib.config

    Args:
        rundate:   Rundate of analysis.
        pipeline:  Pipeline used for analysis.
        session:   Session in analysis.

    Returns:
        Path to configuration file.
    """
    file_vars = config.create_file_vars(rundate, pipeline, session)
    return files.path("config", file_vars=file_vars)


def _dependent_sections(cfg_section, master_cfg=config.where):
    """Find sections that given config section depends on for further configuration

    Looks for keys in the section that has the metadata add_sections set. When such a key is found it looks for
    sections matching the corresponding values.

    Args:
        cfg_section:    Configuration section that will be searched

    Returns:
        Dict:  List of configuration sections
    """
    for key, entry in cfg_section.items():  # Todo: Make ConfigurationSection iterable
        # Check entry for metadata
        if (
            "add_sections" not in entry.meta
            and "add_sections" not in master_cfg.get(key, section=cfg_section.name, default="").meta
        ):
            continue

        # Return sections matching entry values  # Todo: Add .set property to ConfigurationEntry
        sections = set(master_cfg.section_names) & set([s.split(":")[0] for s in entry.list])
        for section in sections:
            yield master_cfg[section]


def _clean_sys_argv():
    """Values in sys.argv that are not valid option values in Where
    """
    reserved_opts = {"id", "session", "profile", "user", "line_profile"}
    return [o for o in sys.argv[1:] if o.startswith("--") and o[2:].split("=")[0] not in reserved_opts]


# Run main function only when running as script
if __name__ == "__main__":
    sys.exit(main())
