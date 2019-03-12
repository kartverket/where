"""Write SISRE buffer file

Description:
------------
By determination of UERE based on Where SISRE results and Terrapos UEE results it is necessary to create a Where buffer
file. In the buffer file the path of the SISRE output file is written after successful processing of a day. The UERE
programs can access this file to check, if something can be processed.
"""
import fcntl

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import log
from where.lib import files


@plugins.register
def sisre_output_buffer(dset):
    """Write SISRE buffer file by appending SISRE output file path

    Args:
        dset:       Dataset, a dataset containing the data.
    """
    with files.open("output_sisre_buffer", file_vars=dset.vars, mode="at") as fid:

        # Allow only one process to hold an exclusive lock for a given file at a given time
        try:
            fcntl.flock(fid, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except IOError:
            log.fatal("flock() failed to hold an exclusive lock.")

        # Append SISRE output file pathes SISRE buffer file
        file_path = files.path(f"output_sisre_2", file_vars=dset.vars)
        fid.write(f"{file_path}\n")

        # Unlock file
        try:
            fcntl.flock(fid, fcntl.LOCK_UN)
        except:
            log.fatal("flock() failed to unlock file.")
