"""Write GNSS dilution of precision results

Description:
------------


"""
# Standard library imports
from collections import namedtuple

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.gnss import gnss
from midgard.writers._writers import get_field, get_header

# Where imports
import where
from where.cleaners.editors.gnss_dop import gnss_dop
from where.lib import config
from where.lib import util

WriterField = namedtuple(
    "WriterField", ["name", "field", "attrs", "dtype", "format", "width", "header", "unit", "description"]
)
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
        description (str):       Description of field
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
    WriterField(
        "date",
        "date",
        (),
        object,
        "%21s",
        19,
        "DATE",
        "YYYY/MM/DD hh:mm:ss",
        "Date in format year, month, day and hour, minute and second",
    ),
    WriterField("mjd", "time", ("gps", "mjd"), float, "%14.6f", 14, "MJD", "", "Modified Julian Day"),
    WriterField("gpsweek", "time", ("gps", "gps_ws", "week"), int, "%5d", 5, "WEEK", "", "GPS week"),
    WriterField(
        "gpssec", "time", ("gps", "gps_ws", "seconds"), float, "%11.3f", 11, "GPSSEC", "second", "GPS seconds"
    ),
    WriterField(
        "sat_avlbl",
        "num_satellite_available",
        (),
        float,
        "%10d",
        10,
        "SAT_AVLBL",
        "",
        "Number of available satellites, which are given as observations",
    ),
    WriterField(
        "sat_used",
        "num_satellite_used",
        (),
        float,
        "%10d",
        10,
        "SAT_USED",
        "",
        "Number of used satellites (healthy, above elevation mask, ...)",
    ),
    WriterField("gdop", "gdop", (), float, "%7.2f", 7, "GDOP", "", "Geometric dilution of precision"),
    WriterField("pdop", "pdop", (), float, "%7.2f", 7, "PDOP", "", "Position dilution of precision"),
    WriterField("tdop", "tdop", (), float, "%7.2f", 7, "TDOP", "", "Time dilution of precision"),
    WriterField("hdop", "hdop", (), float, "%7.2f", 7, "HDOP", "", "Horizontal dilution of precision"),
    WriterField("vdop", "vdop", (), float, "%7.2f", 7, "VDOP", "", "Vertical dilution of precision"),
)


@plugins.register
def gnss_position(dset: "Dataset") -> None:
    """Write GNSS position results


    Args:
        dset:  A dataset containing the data.
    """
    file_path = config.files.path("output_dilution_of_precision", file_vars={**dset.vars, **dset.analysis})

    # Add date and DOPs fields to dataset
    if "date" not in dset.fields:
        dset.add_text(
            "date", 
            val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime], 
            unit="YYYY/MM/DD hh:mm:ss",
            write_level="detail",
        )

    if "pdop" not in dset.fields:
        gnss_dop(dset)

    if "num_satellite_used" not in dset.fields:
        dset.add_float(
            "num_satellite_used",
            val=gnss.get_number_of_satellites(dset.system, dset.satellite, dset.time.gps.datetime),
            write_level="operational",
        )

    # Put together fields in an array as specified by the 'dtype' tuple list
    output_list = list()
    for epoch in dset.unique("time"):
        idx = dset.filter(time=epoch)

        # Append current epoch position solution to final output solution
        output_list.append(tuple([get_field(dset, f.field, f.attrs, f.unit)[idx][0] for f in FIELDS]))

    output_array = np.array(output_list, dtype=[(f.name, f.dtype) for f in FIELDS])

    # Write to disk
    header = get_header(
        FIELDS,
        pgm_version=f"where {where.__version__}",
        run_by=util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else "",
        summary="GNSS dilution of precision results",
    )
    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in FIELDS),
        header=header,
        delimiter="",
        encoding="utf8",
    )
