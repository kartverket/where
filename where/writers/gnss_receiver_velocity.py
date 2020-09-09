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
    # GR-II
    WriterField(
        "sol_val",
        "sol_val",
        (),
        int,
        "%10d",
        10,
        "SOL_VAL",
        "",
        "Solution validation indicator (0=test passed, 1: Not passed)",
    ),
    WriterField(
        "num_satellite_used",
        "num_satellite_used",
        (),
        int,
        "%10d",
        10,
        "SAT_USED",
        "",
        "Number of satellites used in the computation",
    ),
    WriterField("vel_3d", "vel_3d", (), float, "%15.6f", 15, "3DSpeed", "m/s", "Instantaneous 3D speed [m/s]"),
    WriterField("vel_2d", "vel_2d", (), float, "%15.6f", 15, "HSpeed", "m/s", "Instantaneous horizontal speed [m/s]"),
    # GR-III
    ##WriterField("x_vel", "x_vel",  (), float, "%15.6f", 15, "XSpeed" , "m/s", "Instantaneous speed along WGS84 X axis [m/s]"),
    WriterField(
        "y_vel", "y_vel", (), float, "%15.6f", 15, "YSpeed", "m/s", "Instantaneous speed along WGS84 Y axis [m/s]"
    ),
    WriterField(
        "z_vel", "z_vel", (), float, "%15.6f", 15, "ZSpeed", "m/s", "Instantaneous speed along WGS84 Z axis [m/s]"
    ),
    WriterField("e_vel", "e_vel", (), float, "%15.6f", 15, "EastSpeed", "m/s", "Instantaneous east  speed  [m/s]"),
    WriterField("n_vel", "n_vel", (), float, "%15.6f", 15, "NorthSpeed", "m/s", "Instantaneous north  speed  [m/s]"),
    WriterField(
        "u_vel", "u_vel", (), float, "%15.6f", 15, "VerticalSpeed", "m/s", "Instantaneous vertical  speed  [m/s]"
    ),
    # GR-IV
    WriterField("gdop", "gdop", (), float, "%10.2f", 10, "GDOP", "", "Geometric Dilution of Precision"),
    WriterField("pdop", "pdop", (), float, "%10.2f", 10, "PDOP", "", "Position Dilution of Precision"),
    WriterField("hdop", "hdop", (), float, "%10.2f", 10, "HDOP", "", "Horizontal Dilution of Precision"),
    WriterField("vdop", "vdop", (), float, "%10.2f", 10, "VDOP", "", "Vertical Dilution of Precision"),
    ## GR-V
    WriterField(
        "c_xx", "c_xx", (), float, "%10.3f", 10, "COV_XX", "m**2", "Variance of station position X-coordinate"
    ),
    WriterField(
        "c_yy", "c_yy", (), float, "%10.3f", 10, "COV_YY", "m**2", "Variance of station position Y-coordinate"
    ),
    WriterField(
        "c_zz", "c_zz", (), float, "%10.3f", 10, "COV_ZZ", "m**2", "Variance of station position Z-coordinate"
    ),
    WriterField("c_xy", "c_xy", (), float, "%10.3f", 10, "COV_XY", "m**2", "XY-covariance of station position"),
    WriterField("c_xz", "c_xz", (), float, "%10.3f", 10, "COV_XZ", "m**2", "XZ-covariance of station position"),
    WriterField("c_yz", "c_yz", (), float, "%10.3f", 10, "COV_YZ", "m**2", "YZ-covariance of station position"),
)


@plugins.register
def gnss_receiver_velocity(dset: "Dataset") -> None:
    """Write Doppler estimated GNSS receiver velocity results to a text file


    Args:
        dset:  A dataset containing the data.
    """
    file_path = config.files.path("output_receiver_velocity", file_vars=dset.vars)

    # Add date field to dataset
    if "date" not in dset.fields:
        dset.add_text("date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime])

    # Put together fields in an array as specified by the 'dtype' tuple list
    if config.tech.estimate_epochwise.bool:  # Epochwise estimation or over whole time period

        output_list = list()
        for epoch in dset.unique("time"):
            idx = dset.filter(time=epoch)

            # Append current epoch solution to final output solution
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
        summary="GNSS position results",
    )
    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in FIELDS),
        header="\n".join(header),
        delimiter="",
        encoding="utf8",
    )
