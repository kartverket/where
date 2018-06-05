"""Set up a config file for editing data

"""

# Standard library imports
from configparser import ConfigParser
from datetime import datetime
import sys

# Where imports
from midgard.dev import console

from where.lib import config
from where.lib import dependencies
from where.lib import files
from where.lib import log
from where.lib import util


def config_path(rundate, tech):
    """Find path of config file

    Args:
       rundate: The model run date.
       tech:    String, the technique.

    Returns:
        Path: Path to the configuration file for the given date and technique.
    """
    return files.path(
        "config_session",
        file_vars=dict(config.program_vars(rundate=rundate, tech_name=tech), **config.date_vars(date=rundate)),
    )


def add_session_config_file(rundate, tech, overwrite=False):
    """Add an session config file for the current model run

    The session config file is cached for later use.

    Args:
        overwrite:  Overwrite an already existing session file?
    """
    # Check if config file exists
    cfg_path = files.path("config_session")
    if overwrite or not cfg_path.exists():
        # Create and store empty ConfigParser
        session_cfg = ConfigParser()
        log.info("Add session config file for {} at {}", tech.upper(), cfg_path)
        with files.open("config_session", mode="wt") as fid:
            session_cfg.write(fid)

    log.info("Use session config file for {} at {}", tech.upper(), cfg_path)
    dependencies.add(cfg_path)

    # Cache session config for later use
    config.reread()


def add_session_config(session, overwrite=False):
    """Add an session configuration for the given session

    The default contents of the session configuration are copied from the main Where config file. If available, data
    from a section called [{tech}_edit__{user}] (e.g. [vlbi_edit__monitor]) are used. Otherwise the section
    [{tech}_edit] (e.g. [vlbi_edit]) is looked up. If none of these are available, session data are copied from the
    section [default_edit].

    If overwrite is set to False (default) the session section is not changed if it already exists. If overwrite is
    True then the session section will be restored to the default as defined in the main Where config file.

    Args:
        session:    The session for which to add an session file.
        overwrite:  Overwrite an already existing session file?
    """
    # Get session configuration
    if not overwrite and session in config.session.sections:
        return

    # Add session as a config section
    try:
        del config.session[session]
    except KeyError:
        pass
    template = "{}_edit".format(config.analysis.tech.str)
    config.session.update_from_config_section(config.where[template], session)

    lib_cfg_path = files.path("config_session_library")
    lib_cfg = ConfigParser()
    try:
        lib_cfg.read(lib_cfg_path)
        config.session.update_from_dict(lib_cfg[session], section=session, source=str(lib_cfg_path))
        log.info("Using session config {} from library", lib_cfg_path)
    except (FileNotFoundError, KeyError):
        pass

    # Store config to file
    cfg_path = files.path("config_session")
    log.info("Update session config with {} for {} at {}", session, config.analysis.tech.str.upper(), cfg_path)
    with files.open("config_session", mode="wt") as fid:
        print(config.session, file=fid)

    # TODO: Read config from file to be 100% consistent
    # config.reread()   # reread messes up config.analysis


def store_session_config(sessions):
    try:
        fields_to_archive = config.tech.archive_session_fields.list
    except AttributeError:
        fields_to_archive = config.tech.edit_fields_to_archive.list
        log.dev("Config option 'edit_fields_to_archive' is deprecated. Use 'archive_session_fields' instead")
    if not fields_to_archive:
        return

    # Read existing config from library
    lib_cfg = ConfigParser()
    lib_cfg_path = files.path("config_session_library")
    try:
        lib_cfg.read(lib_cfg_path)
    except FileNotFoundError:
        pass

    log.info("Storing session data for sessions {} to {}", ", ".join(sessions), lib_cfg_path)
    for session in sessions:
        if session not in lib_cfg.sections():
            lib_cfg.add_section(session)

        # Add current config settings for the given fields
        for key, entry in config.session[session].items():
            if key in fields_to_archive:
                lib_cfg.set(session, key, entry.str)

        # Store config back to library
        with files.open("config_session_library", mode="wt", create_dirs=True) as fid:
            lib_cfg.write(fid)


def apply_options(rundate, tech, options=None):
    """Apply options to an session configuration file

    The session configuration file is read explicitly by a ConfigParser (instead of using config.get_session_config) to
    avoid saving the metadata that config.get_session_config adds.

    Args:
        rundate (Date):   The model run date.
        tech (String):    The technique.
        options (List):   List of options to apply, by default read from the command line.

    Returns:
        Bool:  True if configuration is changed, False otherwise.

    """
    comment_lines = [
        "Session configuration for {}".format(tech.upper()),
        "Last updated by {} on {}".format(util.get_program_name(), datetime.now().strftime(config.FMT_datetime)),
    ]
    comments = "\n".join(console.fill(c, initial_indent="# ", subsequent_indent="# ") for c in comment_lines) + "\n\n"

    # Use sys.argv as default for options
    if options is None:
        options = sys.argv[1:]

    # Read current configuration file
    cfg_path = config_path(rundate, tech)
    session_cfg = ConfigParser()
    session_cfg.read(cfg_path)

    # Add options
    cfg_changed = ignored_key = False
    for option in [o[10:] for o in options if o.startswith("--session:") and "=" in o]:
        key, value = option.split("=", maxsplit=1)
        *sections, key = key.split(":")
        if not sections:  # Apply to all sections by default
            sections = session_cfg.sections()

        for section in sections:
            if not session_cfg.has_option(section, key):
                log.warn("Section {} of session config file has no key {}. Ignored.", section, key)
                ignored_key = True
                continue
            cfg_changed = cfg_changed or (value != session_cfg[section][key])
            session_cfg.set(section, key, value)

    if cfg_changed:
        with files.open_path(cfg_path, mode="wt") as fid:
            fid.write(comments)
            session_cfg.write(fid)

    return cfg_changed or ignored_key
