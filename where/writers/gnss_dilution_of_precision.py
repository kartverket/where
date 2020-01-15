"""Write GNSS dilution of precision results

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
from where.cleaners.editors.gnss_dop import gnss_dop
from where.lib import config
from where.lib import gnss
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
# # PGM: where 0.21.2/midgard 0.3.0  RUN_BY: NMA  DATE: 20190604 135301 UTC
# #
# #                DATE           MJD WEEK     GPSSEC SAT_AVLBL  SAT_USED   GDOP   PDOP   TDOP   HDOP   VDOP
# # YYYY/MM/DD hh:mm:ss                        second
# # ________________________________________________________________________________________________________
#   2019/02/01 00:00:00  58515.000000 2038 432000.000         8         7   2.28   2.05   1.00   0.96   1.80
#   2019/02/01 00:05:00  58515.003472 2038 432300.000         8         7   2.26   2.04   0.98   0.96   1.80
#   2019/02/01 00:10:00  58515.006944 2038 432600.000         8         7   2.24   2.02   0.96   0.96   1.78
# ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----
#
FIELDS = (
    WriterField("date", "date", (), object, "%21s", 19, "DATE", "YYYY/MM/DD hh:mm:ss"),
    WriterField("mjd", "time", ("gps", "mjd"), float, "%14.6f", 14, "MJD", ""),
    WriterField("gpsweek", "time", ("gps", "gps_ws", "week"), int, "%5d", 5, "WEEK", ""),
    WriterField("gpssec", "time", ("gps", "gps_ws", "seconds"), float, "%11.3f", 11, "GPSSEC", "second"),
    WriterField("sat_avlbl", "num_satellite_available", (), float, "%10d", 10, "SAT_AVLBL", ""),
    WriterField("sat_used", "num_satellite_used", (), float, "%10d", 10, "SAT_USED", ""),
    WriterField("gdop", "gdop", (), float, "%7.2f", 7, "GDOP", ""),
    WriterField("pdop", "pdop", (), float, "%7.2f", 7, "PDOP", ""),
    WriterField("tdop", "tdop", (), float, "%7.2f", 7, "TDOP", ""),
    WriterField("hdop", "hdop", (), float, "%7.2f", 7, "HDOP", ""),
    WriterField("vdop", "vdop", (), float, "%7.2f", 7, "VDOP", ""),
)


@plugins.register
def gnss_position(dset: "Dataset") -> None:
    """Write GNSS position results


    Args:
        dset:  A dataset containing the data.
    """
    file_path = config.files.path("output_dilution_of_precision", file_vars=dset.vars)

    # Add date and DOPs fields to dataset
    if "date" not in dset.fields:
        dset.add_text(
            "date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime], unit="YYYY/MM/DD hh:mm:ss"
        )

    if "pdop" not in dset.fields:
        gnss_dop(dset)

    if "num_satellite_used" not in dset.fields:
        dset.add_float("num_satellite_used", val=gnss.get_number_of_satellites(dset))

    # Put together fields in an array as specified by the 'dtype' tuple list
    output_list = list()
    for epoch in dset.unique("time"):
        idx = dset.filter(time=epoch)

        # Append current epoch position solution to final output solution
        output_list.append(tuple([_get_epoch(dset, idx, f.field, f.attrs, f.unit) for f in FIELDS]))

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


def _get_epoch(dset: "Dataset", idx: np.ndarray, field: "str", attrs: Tuple[str], unit: "str") -> np.ndarray:
    """Get field values of a Dataset specified by the field attributes for a given epoch

    If necessary the unit of the data fields are corrected to the defined 'output' unit.

    Args:
        dset:     Dataset, a dataset containing the data.
        idx:      Dataset indices for selected epoch.
        field:    Field name.
        attrs:    Field attributes (e.g. for Time object: (<scale>, <time format>)).
        unit:     Unit used for output.

    Returns:
        Array with Dataset field values
    """
    f = dset[field]
    for attr in attrs:
        if isinstance(attr, tuple):  # TODO: Workaround. Should be changed e.g. to dset.site_pos.itrs.x with Dataset 3.
            f = f[attr]
            continue

        f = getattr(f, attr)

    # Get first solution of given epoch solutions
    f = f[idx][0]

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
DATE        YYYY/MM/DD hh:mm:ss   Date in format year, month, day and hour, minute and second
MJD                               Modified Julian Day
WEEK                              GPS week
GPSSEC      second                GPS seconds
SAT_AVLBL                         Number of available satellites, which are given as observations
SAT_USED                          Number of used satellites (healthy, above elevation mask, ...)
GDOP                              Geometric dilution of precision
PDOP                              Position dilution of precision
TDOP                              Time dilution of precision
HDOP                              Horizontal dilution of precision
VDOP                              Vertical dilution of precision
"""
    )
    return header
