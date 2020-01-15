"""A parser for reading SP3 orbit files

Example:
--------

    from where import data
    from where import parsers

    # Parse data
    parser = parsers.parse('orbit_sp3', file_path=file_path, rundate=rundate)

    # Create a empty Dataset
    dset = data.Dataset(rundate=rundate, tech=tech, stage=stage, dataset_name=name, dataset_id=0, empty=True)

    # Fill Dataset with parsed data
    parser.write_to_dataset(dset)


Description:
------------
The parser can read orbit data from SP3 file in format c (see :cite:`hilla2010`) or d (see :cite:`hilla2016`).

The file to be parsed should be specified like:

    from where import parsers
    parser = parsers.parse('orbit_sp3', file_path=file_path, rundate=rundate)


"""

# Midgard imports
from midgard.dev import plugins
from midgard.files import dependencies

# Where imports
from where.lib import config
from where.lib import log
from where.parsers import orbit_sp3c, orbit_sp3d


@plugins.register
def get_sp3c_or_sp3d(rundate, file_path=None, **kwargs):
    """Use either OrbitSp3cParser or OrbitSp3dParser for reading orbit files in SP3c or SP3d format

    Firstly the version of SP3 file is read. Based on the read version number it is decided, which Parser should be
    used.

    Args:
        rundate (date):           The model run date.
        file_path (str):          Optional path to orbit-file to parse.
    """
    version = _get_sp3_file_version(file_path)
    dependencies.add(file_path, label="gnss_orbit_sp3")  # MURKS_hjegei: Better solution?

    if version in "ac":
        return orbit_sp3c.OrbitSp3cParser(file_path=file_path, **kwargs)
    elif version.startswith("d"):
        return orbit_sp3d.OrbitSp3dParser(file_path=file_path, **kwargs)
    else:
        log.fatal(f"Unknown SP3 format {version!r} is used in file {file_path}")


def _get_sp3_file_version(file_path):
    """ Get SP3 file version for a given file path

    Args:
        file_path (str):    SP3 file path

    Returns:
        str:            SP3 file version number
    """
    if config.files.empty_file(file_path):
        log.warn(f"File {file_path} is empty.")
    with config.files.open_path(file_path, mode="rt") as fid:
        version = fid.readline().split()[0]

    if len(version) < 2 or version[1] not in "acd":
        log.fatal(f"Unknown SP3 format {version!r} is used in file {file_path}")

    return version[1]
