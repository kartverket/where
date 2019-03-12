"""A parser for reading GNSS RINEX navigation files

Example:
--------

    from where import data
    from where import parsers

    # Parse data
    parser = parsers.parse(file_key=file_key, rundate=rundate, file_vars=file_vars)

    # Create a empty Dataset
    dset = data.Dataset(rundate=rundate, tech=tech, stage=stage, dataset_name=name, dataset_id=0, empty=True)

    # Fill Dataset with parsed data
    parser.write_to_dataset(dset)


Description:
------------

Reads GNSS ephemeris data from RINEX navigation file in format 2.11 (see :cite:`rinex2`) or 3.03 (see
:cite:`rinex3`).


Description:
------------
The parser can read orbit data from SP3 file in format c (see :cite:`hilla2010`) or d (see :cite:`hilla2016`).

The file to be parsed has to be specified for a regular Where-model run with the `file_key`, `rundate` and `file_vars`.
The actual file will then be located in the archive using the files.conf configuration file.


"""

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import gnss
from where.lib import log
from where.parsers import rinex2_nav, rinex3_nav


@plugins.register
def get_rinex2_or_rinex3(rundate, file_vars):
    """Use either Rinex2NavParser or Rinex3NavParser for reading orbit files in format 2.11 or 3.03.

    Firstly the RINEX file version is read. Based on the read version number it is decided, which Parser should be
    used.

    Args:
        rundate (date):           The model run date.
        file_vars (dict):         Dictionary used to replace variables in file name and path.
    """
    # TODO: Due to Geir Arne should handling of file_key be improved.
    version, file_path = gnss.get_rinex_file_version(file_vars["file_key"], file_vars=file_vars)
    if version.startswith("2"):
        return rinex2_nav.Rinex2NavParser(rundate, file_vars["station"])
    elif version.startswith("3"):
        return rinex3_nav.Rinex3NavParser(rundate, file_vars["station"])
    else:
        log.fatal(f"Unknown RINEX format {version} is used in file {file_path}")
