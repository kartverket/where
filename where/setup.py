#!/usr/bin/env python3
""" Utility functions used during start up
"""
# Standard library imports
import editor
from datetime import datetime
import sys

# Midgard imports
from midgard.config import config as mg_config
from midgard.dev import console

# Where imports
from where.lib import config
from where.lib import log
from where.lib import util

_RESERVED_OPTS = {"id", "profile", "user", "line_profile"}


def set_profile(pipeline):
    # Set the correct profile
    profile = util.read_option_value("--profile", default="")
    config.where.profiles = profile.split() + [pipeline]
    config.files.profiles = profile.split() + [pipeline]


def setup_config(rundate, pipeline, *args, **kwargs):
    """Set up configuration for a Where analysis

    """
    set_profile(pipeline)

    # Should a new analysis be started?
    start_new = util.check_options("-N", "--new")

    # Delete an analysis
    if util.check_options("-D", "--delete"):
        from where.tools import delete

        delete.delete_analysis(rundate, pipeline, **kwargs)
        if not start_new:
            raise SystemExit

    # Create configuration of a new analysis
    if start_new or not has_config(rundate, pipeline, *args, **kwargs):
        create_config(rundate, pipeline, *args, **kwargs)
    elif util.check_options("--profile"):  # Warning if --profile option is ignored
        profile_opt = f"--profile={util.read_option_value('--profile', default='')}"
        log.warn(f"Configuration already exists, option '{profile_opt}' ignored")

    # Update configuration based on command line options
    unused_options = update_config(rundate, pipeline, *args, **kwargs)

    # Edit configuration manually
    if util.check_options("-E", "--edit"):
        edit_config(rundate, pipeline, *args, **kwargs)
        unused_options = [opt for opt in unused_options if opt not in ("-E, --edit")]

    # Show current configuration
    if util.check_options("-S", "--show-config"):
        show_config(rundate, pipeline, *args, **kwargs)
        raise SystemExit

    return unused_options


def has_config(rundate, pipeline, *args, **kwargs):
    """Test whether the configuration of a Where analysis exists

    """
    return _config_path(rundate, pipeline, *args, **kwargs).exists()


def create_config(rundate, pipeline, *args, **kwargs):
    """Create the configuration of a Where analysis

    """
    # Create a new configuration and copy all and pipeline sections
    cfg_path = _config_path(rundate, pipeline, *args, **kwargs)
    cfg = mg_config.Configuration(pipeline)
    cfg.update_from_config_section(config.where.all, section=pipeline)
    cfg.update_from_config_section(config.where[pipeline], section=pipeline)
    cfg.write_to_file(cfg_path, metadata=False)

    # Update configuration settings from library
    for section in read_from_library(rundate, pipeline, *args, **kwargs):
        cfg.update_from_config_section(section, section.name)

    # Write updated configuration to file
    cfg.write_to_file(cfg_path, metadata=False)
    log.info(f"Creating new configuration at '{cfg_path}' based on {', '.join(cfg.sources)}")

    # Add new dependent sections from newly created config
    add_sections(rundate, pipeline, *args, **kwargs)

    # Add timestamp and creation note
    add_timestamp(rundate, pipeline, "created", **kwargs)


def update_config(rundate, pipeline, *args, **kwargs):
    """Update the configuration of a Where analysis

    """
    cfg_path = _config_path(rundate, pipeline, *args, **kwargs)
    ts_before = cfg_path.stat().st_mtime

    # Update with command line options
    with mg_config.Configuration.update_on_file(_config_path(rundate, pipeline, *args, **kwargs)) as cfg:
        cfg.master_section = pipeline
        unused_options = cfg.update_from_options(_clean_sys_argv())

    for opt in _RESERVED_OPTS:
        if opt in kwargs:
            unused_options.append(f"--{opt}={kwargs[opt]}")

    # Add timestamp and updated note
    if cfg_path.stat().st_mtime != ts_before:
        add_timestamp(rundate, pipeline, "last update", **kwargs)

    # Add new dependent sections from command line options
    add_sections(rundate, pipeline, *args, **kwargs)

    return unused_options


def edit_config(rundate, pipeline, *args, **kwargs):
    """Update the configuration of a Where analysis

    """
    cfg_path = _config_path(rundate, pipeline, *args, **kwargs)
    ts_before = cfg_path.stat().st_mtime

    # Open config file in an editor
    editor.edit(str(cfg_path))

    if cfg_path.stat().st_mtime != ts_before:
        # Add timestamp and edited note
        add_timestamp(rundate, pipeline, "last update", **kwargs)

        # Add new dependent sections from manual edit
        add_sections(rundate, pipeline, *args, **kwargs)


def add_sections(rundate, pipeline, *args, **kwargs):
    """Update the configuration with sections with settings for models, cleaners etc

    Todo: Figure out how to work with metadata across profiles

    """
    cfg_path = _config_path(rundate, pipeline, *args, **kwargs)
    ts_before = cfg_path.stat().st_mtime

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
    if cfg_path.stat().st_mtime != ts_before:
        add_timestamp(rundate, pipeline, "last update", **kwargs)


def show_config(rundate, pipeline, *args, **kwargs):
    """Show the configuration of a Where analysis

    """
    line = "=" * console.columns()

    # Warn about missing session
    if not has_config(rundate, pipeline, *args, **kwargs):
        log.warn(f"No configuration found for {pipeline.upper()} {rundate.strftime(config.FMT_date)}")

    # Read configuration from file
    else:
        cfg = _read_config(rundate, pipeline, *args, **kwargs)

        # Print configuration to console
        print(line)
        print(f"{pipeline.upper()} {rundate.strftime(config.FMT_date)}\n")
        print(cfg)
        print(f"\nConfig file at {', '.join(cfg.sources)}")


def add_timestamp(rundate, pipeline, timestamp_key, **kwargs):
    """Write or update a timestamp to file

    Args:
        rundate:        Rundate of analysis.
        pipeline:       Pipeline used for analysis.
        session:        Session in analysis.
        timestamp_key:  Key denoting timestamp.
    """
    # Find timestamp file
    file_vars = config.create_file_vars(rundate, pipeline, **kwargs)
    ts_path = config.files.path("timestamp", file_vars=file_vars)

    # Add timestamp with update note to timestamp file
    with mg_config.Configuration.update_on_file(ts_path) as ts_cfg:
        timestamp = f"{datetime.now().strftime(config.FMT_datetime)} by {util.get_program_info()}"
        ts_cfg.update("timestamps", timestamp_key, timestamp, source=__file__)


def read_from_library(rundate, pipeline, *args, **kwargs):
    cfg = _read_config(rundate, pipeline, *args, **kwargs)
    cfg.update_from_options(allow_new=True)
    if not cfg.read_from_library.bool:
        return

    file_vars = config.create_file_vars(rundate, pipeline, **kwargs)
    lib_path = config.files.path("config_library", file_vars=file_vars)
    lib_cfg = mg_config.Configuration.read_from_file("library", lib_path)

    for section in lib_cfg.sections:
        yield section


def store_config_to_library(rundate, pipeline, **kwargs):
    cfg = _read_config(rundate, pipeline, **kwargs)
    if not cfg.write_to_library.bool:
        return

    file_vars = config.create_file_vars(rundate, pipeline, **kwargs)
    lib_path = config.files.path("config_library", file_vars=file_vars)
    lib_cfg = mg_config.Configuration("library")

    for section in cfg.sections:
        for key, entry in section.items():  # Todo: Make ConfigurationSection iterable
            if "library" in entry.meta or "library" in config.where.get(key, section=section.name, default="").meta:
                lib_cfg.update(section.name, key, entry.str, source=entry.source)
                # Todo: Only store entries different from default (issue: profiles?)

    lib_cfg.write_to_file(lib_path)


def _read_config(rundate, pipeline, *args, **kwargs):
    """Read the configuration of a Where analysis from file

    Todo: Add this as a classmethod on Configuration

    Args:
        rundate:   Rundate of analysis.
        pipeline:  Pipeline used for analysis.
        session:   Session in analysis.

    Returns:
        Configuration of Where analysis.
    """
    if not has_config(rundate, pipeline, *args, **kwargs):
        raise FileNotFoundError(f"No configuration found for {pipeline.upper()} {rundate.strftime(config.FMT_date)}")

    cfg = mg_config.Configuration.read_from_file(pipeline, _config_path(rundate, pipeline, *args, **kwargs))
    cfg.master_section = pipeline

    return cfg


def _config_path(rundate, pipeline, *args, **kwargs):
    """The path to the configuration of a Where analysis

    Todo: Move this to lib.config

    Args:
        rundate:   Rundate of analysis.
        pipeline:  Pipeline used for analysis.
        session:   Session in analysis.

    Returns:
        Path to configuration file.
    """
    file_vars = config.create_file_vars(rundate, pipeline, **kwargs)
    return config.files.path("config", file_vars=file_vars)


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
    return [o for o in sys.argv[1:] if o.startswith("--") and o[2:].split("=")[0] not in _RESERVED_OPTS]
