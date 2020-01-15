"""Write GNSS post-fit residual results

Description:
------------


"""
# Standard library imports
from collections import namedtuple
from datetime import datetime
from typing import Tuple

# External library imports
import numpy as np

# Midgard imports
import midgard
from midgard.dev import plugins

# Where imports
import where
from where.lib import config
from where.lib.unit import Unit
from where.lib import util

WriterField = namedtuple("WriterField", ["name", "field", "attrs", "dtype", "format", "width", "header", "unit"])
WriterField.__new__.__defaults__ = (None,) * len(WriterField._fields)
WriterField.__doc__ = """A convenience class for defining a output field for the writer

    Args:
        name  (str):             Unique field name
        field (str):             Dataset field name
        attrs (Tuple[str]):      Field attributes
        dtype (Numpy dtype):     Type of field
        format (str):            Format string
        width (int):             Width of header information
        header (str):            Header information
        unit (str):              Unit of field
    """

# Define fields to plot
#
# # PGM: where 0.21.2/midgard 0.3.0  RUN_BY: NMA  DATE: 20190604 134930 UTC
# #
# #  SAT               EPOCH           MJD WEEK     GPSSEC    AZI   ELEV       RESIDUAL
# #      YYYY/MM/DD hh:mm:ss                        second    deg    deg          meter
# # ___________________________________________________________________________________
#   E07  2018/02/01 00:00:00  58150.000000 1986 345600.000   46.6   22.3         -0.000
#   E11  2018/02/01 00:00:00  58150.000000 1986 345600.000 -107.0    8.0          0.000
#   E12  2018/02/01 00:00:00  58150.000000 1986 345600.000 -114.8   59.0          0.000
# ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+-
#
FIELDS = (
    WriterField("satellite", "satellite", (), object, "%5s", 4, "SAT", ""),
    WriterField("date", "date", (), object, "%21s", 20, "EPOCH", "YYYY/MM/DD hh:mm:ss"),
    WriterField("mjd", "time", ("gps", "mjd"), float, "%14.6f", 14, "MJD", ""),
    WriterField("gpsweek", "time", ("gps", "gps_ws", "week"), int, "%5d", 5, "WEEK", ""),
    WriterField("gpssec", "time", ("gps", "gps_ws", "seconds"), float, "%11.3f", 11, "GPSSEC", "second"),
    WriterField("azimuth", "site_pos", ("azimuth",), float, "%7.1f", 7, "AZI", "deg"),
    WriterField("elevation", "site_pos", ("elevation",), float, "%7.1f", 7, "ELEV", "deg"),
    WriterField("residual", "residual", (), float, "%15.3f", 15, "RESIDUAL", "meter"),
)


@plugins.register
def gnss_residual(dset: "Dataset") -> None:
    """Write GNSS post-fit residual results


    Args:
        dset:  A dataset containing the data.
    """

    file_path = config.files.path("output_residual", file_vars=dset.vars)

    # Add date field to dataset
    if "date" not in dset.fields:
        dset.add_text("date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime])

    # Put together fields in an array as specified by the 'dtype' tuple list
    output_list = list(zip(*(_get_field(dset, f.field, f.attrs, f.unit) for f in FIELDS)))
    output_array = np.array(output_list, dtype=[(f.name, f.dtype) for f in FIELDS])

    # Write to disk
    header = [
        _get_header(dset),
        "".join(f"{f.header:>{f.width}s}" for f in FIELDS),
        "".join(f"{f.unit:>{f.width}s}" for f in FIELDS),
        "_" * sum([f.width for f in FIELDS]),
    ]
    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in FIELDS),
        header="\n".join(header),
        delimiter="",
        encoding="utf8",
    )


def _get_field(dset: "Dataset", field: "str", attrs: Tuple[str], unit: "str") -> np.ndarray:
    """Get field values of a Dataset specified by the field attributes

    If necessary the unit of the data fields are corrected to the defined 'output' unit.

    Args:
        dset:     Dataset, a dataset containing the data.
        field:    Field name.
        attrs:    Field attributes (e.g. for Time object: (<scale>, <time format>)).
        unit:     Unit used for output.

    Returns:
        Array with Dataset field values
    """
    f = dset[field]
    for attr in attrs:
        f = getattr(f, attr)

    # Determine output 'unit'
    # +TODO: Does not work for all fields, because dset.unit() does not except 'time.gps.mjd'.
    if unit.startswith("deg"):
        fieldname = f"{field}.{'.'.join(attrs)}"
        f = f * getattr(Unit, f"{dset.unit(fieldname)[0]}2{unit}")
    # -TODO

    return f


def _get_header(dset: "Dataset") -> str:
    """Get header

    Args:
        dset:  A dataset containing the data.

    Returns:
        Header lines
    """

    pgm = "where " + where.__version__ + "/midgard " + midgard.__version__
    run_by = util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else ""
    file_created = datetime.utcnow().strftime("%Y%m%d %H%M%S") + " UTC"
    header = "PGM: {:s}  RUN_BY: {:s}  DATE: {:s}\n".format(pgm, run_by, file_created)
    header = (
        header
        + """

HEADER      UNIT                  DESCRIPTION
______________________________________________________________________________________________________________________
SAT                               Satellite number
DATE        YYYY/MM/DD hh:mm:ss   Date in format year, month, day and hour, minute and second
MJD                               Modified Julian Day
WEEK                              GPS week
GPSSEC      second                GPS seconds
AZI         degree                Azimuth of satellite in relation to station position
ELEV        degree                Elevation of satellite in relation to station position
RESIDUAL    meter                 Post-fit residuals

"""
    )

    return header
