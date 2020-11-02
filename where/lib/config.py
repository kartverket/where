"""Where library module for handling of Where configuration settings

Example:
--------

    >>> from where.lib import config
    >>> config.where.path.work.path
    PosixPath('/opt/result/where/work')

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
    '{tech}-dataset-{stage}-{yyyy}{mm}{dd}.json'
    >>> config.files.dataset_json.filename.replaced.str
    'vlbi-dataset-{stage}-20091102.json'
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
from functools import lru_cache
import getpass
import pathlib
import sys

# Midgard imports
from midgard.config.config import Configuration
from midgard.config.files import FileConfiguration

# Where imports
from where.lib import enums  # noqa  # Register Where enums


# Base directory of the Where installation
WHERE_DIR = pathlib.Path(__file__).resolve().parent.parent.parent

# Prioritized list of possible names of Where config files
_CONFIG_FILENAMES = dict(
    where=("where_local.conf", "where_pipeline_{pipeline}.conf", "where_pipelines.conf", "where.conf"),
    files=("files_local.conf", "files_pipeline_{pipeline}.conf", "files_pipelines.conf", "files.conf"),
    there=("there_local.conf", "there_pipeline_{pipeline}.conf", "there.conf"),
)

# Prioritized list of possible locations for all Where config files
_CONFIG_DIRECTORIES = (pathlib.Path.cwd(), pathlib.Path.home() / ".where", WHERE_DIR / "config")

# Date format, defined here for consistency
FMT_date = "%Y-%m-%d"

# Datetime format, defined here for consistency
FMT_datetime = "%Y-%m-%d %H:%M:%S"

# Datetime format used in filenames, defined here for consistency
FMT_dt_file = "%Y%m%d-%H%M%S"


def init(rundate, pipeline, **cfg_vars):
    # Update Analysis-configuration
    set_analysis(rundate, pipeline=pipeline, **cfg_vars)

    # Add variables to Files-configuration
    set_file_vars(create_file_vars(rundate, pipeline, **cfg_vars))

    # Read Tech-configuration, and set Where-configuration as parent
    tech.clear()
    cfg_path = files.config.directory.replaced.path / files.config.filename.replaced.path
    tech.update_from_file(cfg_path)
    if pipeline in tech.section_names:
        tech.master_section = pipeline


def read_pipeline(pipeline):
    """Read special pipeline config files"""
    read_where_config(pipeline=pipeline)
    read_files_config(pipeline=pipeline)
    read_there_config(pipeline=pipeline)


def set_analysis(rundate, **cfg_vars):
    """Set analysis configuration
    """
    analysis_vars = program_vars(rundate, **cfg_vars)
    analysis.update_from_dict(analysis_vars, section="config", source="config.set_analysis", allow_new=True)
    analysis.master_section = "config"


def create_file_vars(rundate, pipeline, **cfg_vars):
    """Create variables that can be used with the Files-configuration
    """
    return dict(program_vars(rundate, pipeline, **cfg_vars), **date_vars(rundate))


def set_file_vars(file_vars=None):
    """Add variables to Files-configuration

    Args:
        file_vars (Dict):   Variables that will be made available in files-configuration
    """
    # files.clear_vars()
    files.update_vars({"path_where": str(WHERE_DIR)})
    files.update_vars({"path_{}".format(k): str(v.path) for k, v in where.path.data.items()})
    if file_vars is not None:
        files.update_vars(file_vars)


def program_vars(rundate, pipeline, use_options=True, **other_vars):
    PIPELINE = None if pipeline is None else pipeline.upper()
    prg_vars = dict(rundate="" if rundate is None else rundate.strftime(FMT_date), pipeline=pipeline, PIPELINE=PIPELINE)
    prg_vars.setdefault("user", getpass.getuser().lower())
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

    return prg_vars


@lru_cache()
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
    from where.data.time import Time

    month = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]

    # Create the dict of date variables
    try:
        gpsweek = str(int(Time(date.strftime("%Y-%m-%d %H:%M:%S"), fmt="iso", scale="utc").gps.gps_ws.week))
    except ValueError:
        # gps-scale is not defined for 1970-01-01, which is used by timeseries
        gpsweek = ""
    return dict(
        date=date.strftime("%Y%m%d"),
        yyyy=date.strftime("%Y"),
        ce=date.strftime("%Y")[:2],
        yy=date.strftime("%y"),
        m=str(date.month),
        mm=date.strftime("%m"),
        mmm=month[date.month - 1].lower(),
        MMM=month[date.month - 1].upper(),
        d=str(date.day),
        dd=date.strftime("%d"),
        hh=date.strftime("%H"),
        doy=date.strftime("%j"),
        dow=date.strftime("%w"),
        gpsweek=gpsweek,
    )


def config_paths(cfg_name, **path_vars):
    """Yield all files that contain the given configuration"""
    for file_name in _CONFIG_FILENAMES.get(cfg_name, (f"{cfg_name}.conf",))[::-1]:
        for var, val in path_vars.items():
            file_name = file_name.replace(f"{{{var}}}", val)

        for file_dir in _CONFIG_DIRECTORIES:
            file_path = file_dir / file_name
            if file_path.exists():
                yield file_path
                break


@contextmanager
def update_tech_config(rundate, pipeline, **kwargs):
    file_vars = create_file_vars(rundate, pipeline, **kwargs)
    cfg_path = files.config.directory.replace(**file_vars).path / files.config.filename.replace(**file_vars).path

    with Configuration.update_on_file(cfg_path) as cfg:
        yield cfg

    # Todo: Update timestamp file


def timestamps(rundate, pipeline, **kwargs):
    file_vars = create_file_vars(rundate, pipeline, **kwargs)
    ts_path = files.timestamp.directory.replace(**file_vars).path / files.timestamp.filename.replace(**file_vars).path

    cfg = Configuration.read_from_file("timestamps", ts_path)

    if "timestamps" in cfg.sections:
        return cfg.timestamps.as_dict()
    else:
        return dict()


def read_where_config(**path_vars):
    """Read Where-configuration"""
    where.clear()
    for file_path in config_paths("where", **path_vars):
        where.update_from_file(file_path, interpolate=True)
    where.master_section = "all"


def read_files_config(**path_vars):
    """Read Files-configuration"""
    files.clear()
    for file_path in config_paths("files", **path_vars):
        files.update_from_file(file_path, interpolate=False)
    set_file_vars()


def read_there_config(**path_vars):
    """Read There-configuration"""
    there.clear()
    for file_path in config_paths("there", **path_vars):
        there.update_from_file(file_path, interpolate=True, case_sensitive=True)
    there.master_section = "general"


# Add configurations as module variables
where = Configuration("where")
read_where_config()

files = FileConfiguration("files")
read_files_config()

there = Configuration("there")
read_there_config()

analysis = Configuration("analysis")  # Initialized by config.init later
tech = Configuration("tech")  # Initialized by config.init later
tech.fallback_config = where
