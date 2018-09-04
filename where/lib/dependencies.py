"""Where library module for handling dependencies

Description:
------------

Stores a list of files with a hash/checksum or a timestamp that can be used to detect if a file changes.

Two strategies are available:

- Timestamps: Fast, but not always reliable as timestamps may update without the file actually changing.
- md5 hash/checksum: Slower, since it needs to read through the whole file, but will reliably only trigger when a file
  has changed.
"""

# Standard library imports
import atexit
from datetime import datetime
import json

# Midgard imports
from midgard.config.config import Configuration

# Where imports
from where.lib import files
from where.lib import log

# Variables describing current dependency file
_DEPENDENCY_FILE_VARS = dict()
_CURRENT_DEPENDENCIES = dict()


def init(fast_check=True, **dep_vars):
    """Start a clean list of dependencies

    The dep_vars describe which model run stage the dependency is valid for. These are cached, so after a first
    invocation (as is done in pipelines.run) repeated calls do not need to specify the dep_vars.

    Args:
        fast_check:  Fast check uses timestamps, slow check uses md5 checksums.
        dep_vars:    Variables specifying the model_run_depends-file.
    """
    # Store current dependencies to disk
    write()

    # Update and cache variables
    _DEPENDENCY_FILE_VARS.clear()
    _DEPENDENCY_FILE_VARS["fast_check"] = fast_check
    _DEPENDENCY_FILE_VARS.update(dep_vars)

    # Delete any existing dependency file
    dep_path = files.path("model_run_depends", file_vars=_DEPENDENCY_FILE_VARS)
    try:
        dep_path.unlink()
        log.debug(f"Removing old dependency file {dep_path}")
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
    fast_check = _DEPENDENCY_FILE_VARS["fast_check"]
    for file_path in file_paths:
        file_info = _file_info(file_path, fast_check)
        _CURRENT_DEPENDENCIES[str(file_path)] = file_info
        log.debug(f"Adding dependency: {file_path} ({file_info['checksum']})")


def _file_info(file_path, fast_check):
    """Get file info for a file path

    The contents of the file info depends on whether we are doing a fast check or not.

    Args:
        file_path (String/Path):  File path.

    Returns:
        Dictionary: Info about file.
    """
    file_info = dict(timestamp=files.get_timestamp(file_path))
    if fast_check:
        file_info["checksum"] = file_info["timestamp"]
    else:
        file_info["checksum"] = files.get_md5(file_path)

    return file_info


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

    # Open dependency file or start from a fresh dictionary
    dependency_path = files.path("model_run_depends", file_vars=_DEPENDENCY_FILE_VARS)
    dependencies = Configuration.read_from_file("dependecies", dependency_path)

    # Update dependency information
    for file_path, info in _CURRENT_DEPENDENCIES.items():
        dependencies.update_from_dict(info, section=file_path)

    _CURRENT_DEPENDENCIES.clear()

    # Write to dependency file
    dependencies.write_to_file(dependency_path)


def changed(fast_check=True, **dep_vars):
    """Check if the dependencies of a model run have changed

    Returns True if any of the files in the dependency file have changed, or if the dependency file does not exist.

    Args:
        dep_vars:   Variables specifying the model_run_depends-file.

    Returns:
        Boolean: True if any file has changed or if the dependecy file does not exist.
    """
    # Make sure dependency file exists
    dependency_path = files.path("model_run_depends", file_vars=dep_vars)
    if not dependency_path.exists():
        log.debug(f"Dependency file {dependency_path} does not exist")
        return True

    # Check if any dependencies have changed
    dependencies = Configuration.read_from_file("dependencies", dependency_path)
    for file_path in dependencies.sections:
        previous_checksum = dependencies[file_path].checksum.str
        current_checksum = _file_info(file_path, fast_check=fast_check)["checksum"]
        if current_checksum != previous_checksum:
            log.debug(f"Dependency {file_path} changed from {previous_checksum} to {current_checksum}")
            return True

    return False
