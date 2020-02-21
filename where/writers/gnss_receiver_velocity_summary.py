
"""Write  Doppler GNSS receiver velocity results

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


# =============================================================================================================
# Define fields to plots
#
# Format
# #
# #               EPOCH          MJD WEEK     GPSSEC     satCNT       3DSpeed    HSpeed     XSpeed
# # YYYY/MM/DD hh:mm:ss          Unitless     second     unitless      m/s         m/s        m/s
# # ___________________________________________________________________________________________________________
#  2018/02/01 00:00:00  58150.000000 1986 345600.000    8   0.00256   0.00056    0.0016       ...
#  2018/02/01 00:05:00  58150.003472 1986 345900.000    9   0.00126   0.00026    0.0006       ...
#  2018/02/01 00:10:00  58150.006944 1986 346200.000   13   0.00556   0.00012    0.0006       ...
# ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----
# =============================================================================================================
FIELDS = (
    # GR-I
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
    WriterField("3d_vel", "gnss_spv_vel", (), float, "%10.6f", 10, "3D_VEL", "m/s", ""),
)

# ===========================================================================================
######## enu = (dset.site_pos._itrs2enu @ (dset.site_pos.itrs_pos - ref_pos)[:,:,None])[:,:,0]
# ===========================================================================================
@plugins.register
def gnss_receiver_velocity_summary(dset: "Dataset") -> None:
    """Write Doppler estimated GNSS receiver velocity results to a text file


    Args:
        dset:  A dataset containing the data.
    """
    # File name generation
    file_path = config.files.path("output_receiver_velocity_summary", file_vars=dset.vars)

    # Handle GR-I: date field to dataset
    if "date" not in dset.fields:
        dset.add_text("date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime])

    # Put together fields in an array as specified by the 'dtype' tuple list
    if config.tech.estimate_epochwise.bool:  # Epochwise estimation or over whole time period

        output_list = list()
        for epoch in dset.unique("time"):
            idx = dset.filter(time=epoch)

            # Append current epoch position solution to final output solution
            output_list.append(tuple([get_field(dset, f.field, f.attrs, f.unit)[idx][0] for f in FIELDS]))

    else:
        # Get position solution for first observation
        idx = np.squeeze(np.array(np.nonzero(dset.time.gps.mjd)) == 0)  # first observation -> TODO: Better solution?
        output_list = [tuple([get_field(dset, f.field, f.attrs, f.unit)[idx][0] for f in FIELDS])]

    # define the output atrray
    output_array = np.array(output_list, dtype=[(f.name, f.dtype) for f in FIELDS])

    # Write to disk
    header = get_header(
        FIELDS,
        pgm_version=f"where {where.__version__}",
        run_by=util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else "",
        summary="GNSS station velocity results",
    )

    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in FIELDS),
        header="\n".join(header),
        delimiter="",
        encoding="utf8",
    )
