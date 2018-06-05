#!/usr/bin/env python3
"""Set-up configuration files before running a Where analysis

Script:
-------

The script can be run independently of Where as follows::

    ./setup_config.py date technique [options]

with

================= =============================================================
Option            Description
================= =============================================================
date              The model run date in the format ``<year month day>``.
technique         One of the recognized techniques, e.g. vlbi, slr, gnss.

-h, --help        Show this help message and exit.
--debug, ...      Show additional debug information. Other flags such as
                  --all, --debug, --info, --warn, --error, --fatal, --none
                  are also allowed, and will show differing amounts of
                  information as the program runs.
--version         Show version information and exit.
================= =============================================================

By default, the configuration is taken from the general :file:`where.conf`-file.
Individual options can be overridden by specifying ``--option=value`` on the
command line. It is simpler and more common to run the main Where program, which
handles all necessary dependencies (see :doc:`user_guide_where`).


Description:
------------

The script sets up the configuration file for a given technique and model run
date. If no configuration file already exists, a default configuration is taken
from the general :file:`where.conf`-file. If the configuration file exists
beforehand, only individual options specified on the command line are applied.
"""

# Standard library imports
from configparser import ConfigParser
from datetime import datetime
import sys

# Where imports
from midgard.dev import console

import where
from where.lib import config
from where.lib import files
from where.lib import log
from where.lib import util


def main():
    """Parse command line arguments and set up a technique specific configuration file
    """

    # Initialize
    rundate, tech = util.parse_args("date", "string")

    # Run program
    if not is_configured(rundate, tech):
        default_config(rundate, tech)
    apply_options(rundate, tech)


def config_path(rundate, tech):
    """Find path of config file

    Args:
       rundate: The model run date.
       tech:    String, the technique.

    Returns:
        Path: Path to the configuration file for the given date and technique.
    """
    return files.path(
        "config",
        file_vars=dict(config.program_vars(rundate=rundate, tech_name=tech), **config.date_vars(date=rundate)),
    )


def is_configured(rundate, tech):
    """Checks if a configuration file exists for the given run date and tech

    Args:
       rundate: The model run date.
       tech:    String, the technique.

    Returns:
        Bool: True if the configuration file exists, False otherwise.
    """
    return config_path(rundate, tech).exists()


def default_config(rundate, tech, write_to_file=True):
    """Create a default configuration file based on where.conf

    Args:
       rundate: The model run date.
       tech:    String, the technique.
    """
    comments = _format_comments(
        "Configuration for {} based on {}.".format(tech.upper(), ", ".join(config.where.sources)),
        "Written by {} on {}\n".format(util.get_program_name(), datetime.now().strftime(config.FMT_datetime)),
    )

    # Figure out which sections of the main configuration to copy.
    profile = util.read_option_value("--profile", default="")
    config.where.profiles = [profile, tech] if profile else [tech]
    sections = [tech] + _check_for_trigger_words(tech)

    # Copy configuration from the appropriate section in where.conf
    default_cfg = ConfigParser()
    for section in sections:
        default_cfg.add_section(section)
        cfg_section = config.where[section] if section in config.where.sections else {}
        for key, value in cfg_section.items():
            if "__" not in key:  # Do not include profile specific configuration
                default_cfg.set(section, key, value.str)

        # Special profile configuration
        if profile:
            for key, value in cfg_section.items():
                if key.startswith(profile + "__"):
                    default_cfg.set(section, key[len(profile) + 2:], value.str)

    # Common and dynamic options added to the tech part of the configuration
    cfg_all = config.where["all"] if "all" in config.where.sections else {}
    for key, value in cfg_all.items():
        if key not in default_cfg[tech]:
            default_cfg.set(tech, key, value.str)

    default_cfg.set(tech, "timestamp", datetime.now().strftime(config.FMT_dt_file))
    default_cfg.set(tech, "version", where.__version__)

    # Write configuration file
    if write_to_file:
        with files.open_path(config_path(rundate, tech), create_dirs=True, mode="wt") as fid:
            fid.write(comments)
            default_cfg.write(fid)

    return default_cfg


def update_config(rundate, tech):
    """Update options in a configuration file, based on the master configuration

    Args:
        rundate (Date):  The model run date.
        tech (String):   The technique.

    Returns:
        Boolean: True if configuration changed, False otherwise.
    """
    cfg_changed = False
    comments = _format_comments(
        "Configuration for {}".format(tech.upper()),
        "Last updated by {} on {}".format(util.get_program_name(), datetime.now().strftime(config.FMT_datetime)),
    )

    # Read default configuration
    default_cfg = default_config(rundate, tech, write_to_file=False)

    # Read current configuration file
    cfg_path = config_path(rundate, tech)
    tech_cfg = ConfigParser()
    tech_cfg.read(cfg_path)

    # Copy missing configurations
    for section in default_cfg.sections():
        if section in tech_cfg:
            for key, value in default_cfg[section].items():
                if key not in tech_cfg[section]:
                    cfg_changed = True
                    tech_cfg.set(section, key, value)
        else:
            cfg_changed = True
            tech_cfg.add_section(section)
            for key, value in default_cfg[section].items():
                tech_cfg.set(section, key, value)

    # Write configuration file
    if cfg_changed:
        with files.open_path(cfg_path, mode="wt") as fid:
            fid.write(comments)
            tech_cfg.write(fid)

    return cfg_changed


def apply_options(rundate, tech, options=None):
    """Apply additional options to a configuration file

    Args:
        rundate: The model run date.
        tech:    String, the technique.
        options: List of options to apply, by default read from command line.

    Returns:
        Bool: True if configuration changed, False otherwise.
    """
    reserved_keys = {
        "profile", "id", "only_session", "show_profile", "profile_output", "line_profile", "user"
    }  # TODO: better solution?

    comments = _format_comments(
        "Configuration for {}".format(tech.upper()),
        "Last updated by {} on {}".format(util.get_program_name(), datetime.now().strftime(config.FMT_datetime)),
    )

    # Use sys.argv as default for options
    if options is None:
        options = sys.argv[1:]

    # Read current configuration file
    cfg_path = config_path(rundate, tech)
    tech_cfg = ConfigParser()
    tech_cfg.read(cfg_path)

    # Add options
    cfg_changed = ignored_key = False
    for option in [o[2:] for o in options if o.startswith("--") and "=" in o]:
        if option.startswith("session:"):
            continue

        section = tech
        key, value = option.split("=", maxsplit=1)
        if ":" in key:
            section, key = key.split(":")
        if key in reserved_keys:
            continue
        if not tech_cfg.has_option(section, key):
            log.warn("{} has no key {}. Ignored", section, key)
            ignored_key = True
            continue
        cfg_changed = cfg_changed or (value != tech_cfg[section][key])
        tech_cfg.set(section, key, value)

    # Write configuration file
    if cfg_changed:
        with files.open_path(cfg_path, mode="wt") as fid:
            fid.write(comments)
            tech_cfg.write(fid)

    return cfg_changed or ignored_key


def _format_comments(*comments):
    """Format comments with #s and proper linebreaks

    One or more comment strings are taken as input. A line break is guaranteed between different strings, but may also
    be inserted in the middle of a given string if necessary.

    Args:
        comments:  Strings with comments.

    Returns:
        String: Formatted comments.
    """
    return "\n".join(console.fill(c, initial_indent="# ", subsequent_indent="# ") for c in comments) + "\n\n"


def _check_for_trigger_words(*sections_to_check):
    """Finds sections in config file specifying configurations for individual models etc

    Looks at keys that contains trigger words like ``models`` or ``estimate`` and checks if there are sections
    available matching the values of those keys.

    Args:
        sections_to_check:   Strings, sections containing keys with trigger words.

    Returns:
        List: List of strings with names of additional sections in config file.
    """
    trigger_words = {"models", "estimate", "output", "process"}
    sections = []

    for section in sections_to_check:
        for key, value in config.where[section].items():
            if "__" in key:
                continue
            words_in_key = set(w.lower() for w in key.replace("_", " ").split())
            if words_in_key & trigger_words:
                sections += [v.split(":")[0] for v in value.list]

    return sorted(set(sections) & set(config.where.sections))


# Run main function only when running as script
if __name__ == "__main__":
    sys.exit(main())
