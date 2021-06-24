"""Write GNSS post-fit residual results

Description:
------------


"""
# Standard library imports
from collections import namedtuple

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.writers._writers import get_field, get_header

# Where imports
import where
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
    WriterField("satellite", "satellite", (), object, "%5s", 4, "SAT", "", "Satellite number"),
    WriterField(
        "date",
        "date",
        (),
        object,
        "%21s",
        20,
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
        "azimuth",
        "site_pos",
        ("azimuth",),
        float,
        "%7.1f",
        7,
        "AZI",
        "deg",
        "Azimuth of satellite in relation to station position",
    ),
    WriterField(
        "elevation",
        "site_pos",
        ("elevation",),
        float,
        "%7.1f",
        7,
        "ELEV",
        "deg",
        "Elevation of satellite in relation to station position",
    ),
)


@plugins.register
def gnss_residual(dset: "Dataset") -> None:
    """Write GNSS post-fit residual results

    Args:
        dset:  A dataset containing the data.
    """

    file_path = config.files.path("output_residual", file_vars={**dset.vars, **dset.analysis})
    
    # Update WriterField depending on used pipeline
    fields_def = list(FIELDS)
    fields_def.append(WriterField(
                        "residual", 
                        "residual", 
                        (), 
                        float, 
                        "%15.4f" if dset.vars["pipeline"] == "gnss_vel" else "%15.3f", 
                        15, 
                        "RESIDUAL", 
                        "meter/second" if dset.vars["pipeline"] == "gnss_vel" else "meter", 
                        "Post-fit residual",
                      )
    )

    # Add date field to dataset
    if "date" not in dset.fields:
        dset.add_text("date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime], write_level="detail")

    # Put together fields in an array as specified by the 'dtype' tuple list
    output_list = list(zip(*(get_field(dset, f.field, f.attrs, f.unit) for f in fields_def)))
    output_array = np.array(output_list, dtype=[(f.name, f.dtype) for f in fields_def])

    # Write to disk
    header = get_header(
        fields_def,
        pgm_version=f"where {where.__version__}",
        run_by=util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else "",
        summary="GNSS post-fit residual results",
    )

    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in fields_def),
        header=header,
        delimiter="",
        encoding="utf8",
    )
