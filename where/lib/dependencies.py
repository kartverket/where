"""Where library module for handling dependencies

Example:
--------

from where.lib import files
with files.open('eopc04_iau', mode='rt') as fid:
    for line in fid:
        print(line.strip())

Description:
------------

This module handles opening of files. All regular Where files should be registered in the Where file list, and then
opened using files.open. Other files could be opened using files.open_path, although this should be the exception,
rather than the rule.

The files.open and files.open_path functions are both wrappers around the built-in open function, and behave mainly
similar. In particular, they accept all the same keyword arguments (like for instance mode). Furthermore, to make sure
files are properly closed they should normally be used with a context manager as in the example above.


$Revision: 15267 $
$Date: 2018-06-06 01:18:55 +0200 (Wed, 06 Jun 2018) $
$LastChangedBy: hjegei $

"""

# Standard library imports
import atexit
from datetime import datetime
import json
import pathlib

# Where imports
from where.lib import files
from where.lib import log

# Variables describing current dependency file
_DEPENDENCY_FILE_VARS = dict()
_CURRENT_DEPENDENCIES = dict()


def init(**dep_vars):
    """Start a clean list of dependencies

    The dep_vars describe which model run stage the dependency is valid for. These are cached, so after a first
    invocation (as is done in pipelines.run) repeated calls do not need to specify the dep_vars.

    Args:
        dep_vars:   Variables specifying the model_run_depends-file.
    """
    # Store current dependencies to disk
    write()

    # Update and cache variables
    _DEPENDENCY_FILE_VARS.clear()
    _DEPENDENCY_FILE_VARS.update(dep_vars)

    # Delete any existing dependency file
    dep_path = files.path("model_run_depends", file_vars=_DEPENDENCY_FILE_VARS)
    try:
        dep_path.unlink()
        log.debug("Removing old dependency file {}", dep_path)
    except FileNotFoundError:
        pass  # If dependency file does not exist, we do not do anything

    # Register _write in case program exits without writing all dependencies to disk
    atexit.register(_write)


def add(*file_paths):
    """Add a list of files to the list of dependencies

    Records the current time stamp of the files specified by file paths, and stores as dependencies on the dependency
    file.

    Before adding dependencies, a call to `init_dependencies` has to be done, to set up where to store the
    dependencies.

    Args:
        file_paths:   List of file paths.
    """
    # Ignore dependency if no dependency variables are available (init_dependecies has not been called)
    if not _DEPENDENCY_FILE_VARS:
        return

    # Add or update dependency information
    for file_path in file_paths:
        timestamp = files.get_timestamp(file_path)
        _CURRENT_DEPENDENCIES[str(file_path)] = timestamp
        log.debug(f"Adding dependency: {file_path} modified at {timestamp}")


def write():
    """Write dependencies to file
    """
    atexit.unregister(_write)
    _write(write_as_crash=False)


def _write(write_as_crash=True):
    """Write dependencies to file

    This function is called either when starting a new list of dependencies (with a call to `init`) or when the program
    exits (including with an error). If `write_as_crash` is True, a special dependency is stored that will force
    `changed` to return True. This will in particular make sure that a stage is rerun if it crashed the previous time
    it ran.

    Args:
        write_as_crash (Boolean):   Whether to note that the current dependendee crashed.
    """
    # Ignore dependency if no dependency variables are available (init_dependecies has not been called)
    if not _DEPENDENCY_FILE_VARS:
        return

    # Store timestamp of crash, this will also force the current stage to be rerun next time
    if write_as_crash:
        _CURRENT_DEPENDENCIES["__crashed__"] = datetime.now().isoformat()

    # No need to open and close files if there are no dependencies to store
    if not _CURRENT_DEPENDENCIES:
        return

    # Open dependency file
    try:
        with files.open("model_run_depends", file_vars=_DEPENDENCY_FILE_VARS, mode="rt", write_log=False) as fid:
            dependencies = json.load(fid)
    except FileNotFoundError:
        dependencies = dict()

    # Update dependency information
    dependencies.update(_CURRENT_DEPENDENCIES)
    _CURRENT_DEPENDENCIES.clear()

    # Write to dependency file
    try:
        with files.open("model_run_depends", file_vars=_DEPENDENCY_FILE_VARS, mode="wt") as fid:
            json.dump(dependencies, fid)
    except FileNotFoundError:
        log.warn(
            "Not able to write dependencies to {}", files.path("model_run_depends", file_vars=_DEPENDENCY_FILE_VARS)
        )


def changed(**dep_vars):
    """Check if the dependencies of a model run have changed

    Returns True if any of the files in the dependency file have changed, or if the dependency file does not exist.

    Args:
        dep_vars:   Variables specifying the model_run_depends-file.

    Returns:
        Boolean: True if any file has changed or if the dependecy file does not exist.
    """
    # Open dependency file
    try:
        with files.open("model_run_depends", file_vars=dep_vars, mode="rt") as fid:
            dependencies = json.load(fid)
    except FileNotFoundError:
        log.debug("Dependency file {} does not exist", files.path("model_run_depends", file_vars=dep_vars))
        return True

    # Check if any dependencies have changed
    for file_path, timestamp in dependencies.items():
        if timestamp != files.get_timestamp(file_path):
            log.debug("Dependency {} changed from {} to {}", file_path, timestamp, files.get_timestamp(file_path))
            return True

    return False
