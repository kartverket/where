"""Where library module for handling of Where configuration settings

Example:
--------

    >>> from where.lib import config
    >>> config.where.path.work.path
    PosixPath('/opt/result/where/work')
    >>> config.analysis.rundate
    ConfigurationEntry('rundate', '2009-11-02')
    >>> config.analysis.rundate.str
    '2009-11-02'
    >>> config.analysis.rundate.date
    datetime.date(2009, 11, 2)
    >>> config.session.pos_models
    ConfigurationEntry('pos_models', 'ocean_tides, solid_tides, solid_pole_tides')
    >>> config.session.pos_models.tuple
    ('ocean_tides', 'solid_tides', 'solid_pole_tides')


Description:
------------

This module is used to read Where configuration settings. We first try to read configuration settings from the current
working directory, then from Where's config directory (see `_config_directories()`).  The main configuration files
should be called where.conf or where_<pipeline>.conf. Personal changes to the config can be done in a file called
where_local.conf (see `_config_filenames()`).

The configuration is split into four parts:

- `files`    - Information about the files Where knows about. Read from `files.conf`.
- `where`    - Main configuration of Where. Read from `where.conf` and `where_local.conf`.
- `analysis` - Information about the running program. Mainly based on command line input.
- `session`  - Configuration of the current analysis. Copied from the `where` configuration and stored in the work dir.

Each of these are further split into sections, and each section consists of `key=value`-pairs. To read a configuration
entry, use `config.part.section.key`, for instance `config.where.path.work` reads the key `work` in the `path`-section
of the `where` configuration. To actually use a configuration entry you should convert it to the required data type
using one of the properties `str`, `int`, `float`, `bool`, `list`, `tuple`, `dict`, `date`, `datetime` or `path`.

The `files` configuration also has some variables defined based on the current analysis (if `config.init` has been
run). These can be used to replace the values of the configuration entries. For instance,

    >>> config.files.dataset_json.filename.str
    '{$tech}-dataset-{$stage}-{$yyyy}{$mm}{$dd}.json'
    >>> config.files.dataset_json.filename.replaced.str
    'vlbi-dataset-{$stage}-20091102.json'
    >>> config.files.dataset_json.filename.replace(stage='read').str
    'vlbi-dataset-read-20091102.json'

As the last example shows, it is also possible to define the value of some of these variables dynamically, by using the
function `replace` instead of the property `replaced`. In either case, a new configuration entry is returned
which can be converted in the usual way (above we convert to `string`).

Each configuration entry has an associated `source` indicating the origin of the value. This may be a filename, a
descriptive string and/or indications of how the value has been modified. Some examples are:

    >>> # Value from file
    >>> config.files.dataset_json.filename.source
    '/home/where/config/files.conf'
    >>> # Value from file, replaced with values in parenthesis
    >>> config.files.dataset_json.filename.replace(stage='read').source
    '/home/hjegei/where/config/files.conf (tech=vlbi,stage=read,yyyy=2009,mm=11,dd=02)'
    >>> # Value from command line
    >>> config.analysis.id.source
    'command line'

"""

# Standard library imports
from contextlib import contextmanager
import getpass
import pathlib
import sys

# Where imports
from midgard.config import config as mg_config
from where.lib import enums  # noqa  # Register Where enums


# Base directory of the Where installation
WHERE_DIR = pathlib.Path(__file__).parent.parent.parent.resolve()

# Possible names of Where config files
_CONFIG_FILENAMES = dict(
    where=("where_local.conf", "where_pipeline_*.conf", "where_pipelines.conf", "where.conf"),
    files=("files_local.conf", "files_pipeline_*.conf", "files_pipelines.conf", "files.conf"),
    there=("there_local.conf", "there.conf"),
)

# Possible locations for all Where config files
_CONFIG_DIRECTORIES = (pathlib.Path.cwd(), WHERE_DIR / "config")

# Date format, defined here for consistency
FMT_date = "%Y-%m-%d"

# Datetime format, defined here for consistency
FMT_datetime = "%Y-%m-%d %H:%M:%S"

# Datetime format used in filenames, defined here for consistency
FMT_dt_file = "%Y%m%d-%H%M%S"


def init(rundate, tech_name, session, **cfg_vars):
    # Update Analysis-configuration
    set_analysis(rundate, tech=tech_name, session=session, **cfg_vars)

    # Temporarily keep analysis variables in program config
    prg_vars = program_vars(rundate, tech_name, session=session, **cfg_vars)
    program.update_from_dict(prg_vars, section="program", source="config.init", allow_new=True)
    program.master_section = "program"

    # Add variables to Files-configuration
    set_file_vars(create_file_vars(rundate, tech_name, session=session, **cfg_vars))

    # Read Tech-configuration, and set Where-configuration as parent
    # TODO: Rename to analysis?
    tech.clear()
    cfg_path = files.config.directory.replaced.path / files.config.filename.replaced.path
    tech.update_from_file(cfg_path)
    if tech_name in tech.sections:
        tech.master_section = tech_name


def set_analysis(rundate, **cfg_vars):
    """Set analysis configuration

    TODO; The analysis config should eventually replace the program config
    """
    analysis_vars = program_vars(rundate, cfg_vars["tech"], **cfg_vars)
    analysis.update_from_dict(analysis_vars, section="config", source="config.set_analysis", allow_new=True)
    analysis.master_section = "config"


def create_file_vars(rundate, pipeline, session, **cfg_vars):
    """Create variables that can be used with the Files-configuration
    """
    return dict(program_vars(rundate, pipeline, session=session, **cfg_vars), **date_vars(rundate))


def set_file_vars(file_vars=None):
    """Add variables to Files-configuration

    Args:
        file_vars (Dict):   Variables that will be made available in files-configuration
    """
    files.clear_vars()
    files.update_vars({"path_where": str(WHERE_DIR)})
    files.update_vars({"path_{}".format(k): str(v.path) for k, v in where.path.data.items()})
    if file_vars is not None:
        files.update_vars(file_vars)


def reread():
    # TODO: Can this function be deleted?
    rundate = program.rundate.date
    tech_name = program.tech.str
    session = program.session.str
    init(rundate=rundate, tech_name=tech_name, session=session)


def program_vars(rundate, tech_name, session, use_options=True, **other_vars):
    prg_vars = dict(rundate="" if rundate is None else rundate.strftime(FMT_date), tech=tech_name, session=session)
    prg_vars.setdefault("user", getpass.getuser())
    prg_vars.setdefault("id", "")

    # Add other variables
    prg_vars.update(other_vars)

    # Update variables from the command line, only those already in prg_vars
    if use_options:
        options = sys.argv[1:]
        for option in options:
            if not (option.startswith("--") and "=" in option):
                continue

            opt_key, _, opt_value = option[2:].partition("=")
            if opt_key in prg_vars:
                prg_vars[opt_key] = opt_value

    # Update id to make sure it starts with an underscore, `_`
    if prg_vars["id"] and not prg_vars["id"].startswith("_"):
        prg_vars["id"] = "_{}".format(prg_vars["id"])

    return prg_vars


def date_vars(date):
    """Construct a dict of date variables

    From a given date, construct a dict containing all relevant date variables. This dict can be used to for instance
    replace variables in file names.

    Examples:
        >>> from datetime import date
        >>> date_vars = date_vars(date(2009, 11, 2))
        >>> sorted(date_vars.items())    # doctest: +NORMALIZE_WHITESPACE
        [('MMM', 'NOV'), ('ce', '20'), ('d', '2'), ('dd', '02'), ('dow', '1'), ('doy', '306'), ('gpsweek', '1556'),
         ('m', '11'), ('mm', '11'), ('mmm', 'nov'), ('yy', '09'), ('yyyy', '2009')]

    Args:
        date (Date/Datetime):      The date.

    Returns:
        Dict: Dictionary with date variables for the given date.
    """
    if date is None:
        return dict()

    # Import Time locally to avoid circular imports
    from where.lib.time import Time

    # Create the dict of date variables
    return dict(
        yyyy=date.strftime("%Y"),
        ce=date.strftime("%Y")[:2],
        yy=date.strftime("%y"),
        m=str(date.month),
        mm=date.strftime("%m"),
        mmm=date.strftime("%b").lower(),
        MMM=date.strftime("%b").upper(),
        d=str(date.day),
        dd=date.strftime("%d"),
        hh=date.strftime("%H"),
        doy=date.strftime("%j"),
        dow=date.strftime("%w"),
        gpsweek=str(int(Time(date.strftime("%Y-%m-%d %H:%M:%S"), format="iso").gpsweek)),
    )


def _config_paths(cfg_name):
    for file_name in _CONFIG_FILENAMES[cfg_name][::-1]:
        for file_dir in _CONFIG_DIRECTORIES:
            file_paths = sorted(pathlib.Path(file_dir).glob(file_name))
            if file_paths:
                for file_path in file_paths:
                    yield file_path
                break


@contextmanager
def update_tech_config(rundate, tech, session, **kwargs):
    file_vars = create_file_vars(rundate, tech, session, **kwargs)
    cfg_path = files.config.directory.replace(**file_vars).path / files.config.filename.replace(**file_vars).path

    with mg_config.Configuration.update_on_file(cfg_path) as cfg:
        yield cfg

    # Todo: Update timestamp file


def reset_config():
    # Where-configuration
    where.clear()
    for file_path in _config_paths("where"):
        where.update_from_file(file_path, interpolate=True)
    where.master_section = "all"

    # Files-configuration
    files.clear()
    for file_path in _config_paths("files"):
        files.update_from_file(file_path, interpolate=False)
    set_file_vars()

    # There-configuration
    there.clear()
    for file_path in _config_paths("there"):
        there.update_from_file(file_path, interpolate=True)
    there.master_section = "general"

    # Analysis, Tech and Session-configurations are initialised by config.init later
    analysis.clear()
    tech.clear()
    tech.parent_config = where
    session.clear()


# Add configurations as module variables
where = mg_config.Configuration("where")
files = mg_config.Configuration("files")
program = mg_config.Configuration("program")
analysis = mg_config.Configuration("analysis")
tech = mg_config.Configuration("tech")
session = mg_config.Configuration("session")
there = mg_config.Configuration("there")
reset_config()
