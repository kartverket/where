"""Delete a Where analysis

Usage:

    {exe:tools} delete <date> <pipeline> [--session=<session>] [options]

The tool requires a date given in the format `<year month day>` (for example
2015 8 4).

In addition, one pipeline must be specified. See below for available pipelines.

===================  ===========================================================
Pipeline             Description
===================  ===========================================================
{pipelines_doc:Delete}
===================  ===========================================================


Description:
------------

The delete tool deletes all files for a given model run from the work
directory.


Examples:
---------

Delete the VLBI analysis for August 4 2015:

    {exe:tools} delete 2015 8 4 -v --session=XA


Current Maintainers:
--------------------

{maintainers}

Version: {version}

"""

# Standard library imports
import pathlib
import shutil

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import files
from where.lib import log


@plugins.register
def delete_analysis(rundate: "date", pipeline: "pipeline", session: "option" = ""):  # typing: ignore
    """Delete working directory for a given model run date

    Args:
        rundate: The model run date.
    """
    file_vars = config.create_file_vars(rundate, pipeline, session=session)
    work_directory = files.path("directory_work", file_vars=file_vars)
    log.info(f"Deleting '{work_directory}'")
    _warn_about_cwd_deleted(work_directory)
    try:
        shutil.rmtree(work_directory)
    except FileNotFoundError:
        log.warn(f"'{work_directory}' does not exist. Nothing to delete")


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
