#!/usr/bin/env python3
"""Archive a Where analysis

Description:
------------

The :func:`archive_analysis` function moves all files for a given model run from
the work directory to the archive. Locations of the work and archive directories
are defined in the :file:`files.conf` configuration file.
"""

# Standard library imports
import pathlib
import shutil

# Midgard imports
from midgard.config import config as mg_config

# Where imports
from where.lib import config
from where.lib import files
from where.lib import log


def archive_analysis(rundate, pipeline, session=None):
    """Move an analysis to the archive for a given model run date

    The archive directory is timestamped. If there are several techniques with timestamps, the earliest one is used for
    the archive directory.

    Args:
        rundate: The model run date.
        tech:    Technique

    Returns:
        String: The timestamp of the archived analysis.
    """
    file_vars = config.create_file_vars(rundate, pipeline, session)
    work_directory = files.path("directory_work", file_vars=file_vars)
    ts_path = files.path("timestamp", file_vars=file_vars)
    try:
        cfg = mg_config.Configuration.read_from_file("timestamps", ts_path)
    except FileNotFoundError:
        log.warn(f"No analysis found at '{work_directory}'. Nothing moved to archive")
        return

    # Move the analysis from the work path to the archive
    timestamp = cfg.timestamps.created.tuple[0]
    archive_directory = files.path("directory_archive", file_vars=dict(file_vars, timestamp=timestamp))
    log.info(f"Moving '{work_directory}' to '{archive_directory}'")
    _warn_about_cwd_deleted(work_directory)
    shutil.move(work_directory, archive_directory)


def _warn_about_cwd_deleted(directory):
    """Warn about the current working directory being deleted

    Deleting the current working directory can lead to weird bugs when continuing to work in the
    terminal. Unfortunately, we cannot change the working directory of the terminal from within Where, as that is a
    parent process. For now, we simply print a warning, and instructions about changing to an existing directory.

    Args:
        directory (Path):  Directory that will be deleted
    """
    cwd = pathlib.Path.cwd()
    directory = directory.resolve()
    if directory == cwd or directory in cwd.parents:
        log.warn(f"Current working directory '{cwd}' is being deleted")
        log.warn(f"Do 'cd {directory.parent}' to keep working")
