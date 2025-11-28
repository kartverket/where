"""Write GNSS velocity results

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
from where.postprocessors.gnss_velocity_fields import gnss_velocity_fields

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
# 
#                DATE           MJD WEEK     GPSSEC     ECEF [X]     ECEF [Y]     ECEF [Z]    SIGMA_X    SIGMA_Y 
# YYYY/MM/DD hh:mm:ss                        second          m/s          m/s          m/s        m/s        m/s 
# _______________________________________________________________________________________________________________
# 2019/02/01 00:00:00  58515.000000 2038 432000.000      -0.0001       0.0054       0.0060     0.0055     0.0030 
# 2019/02/01 00:20:00  58515.013889 2038 433200.000      -0.0146       0.0069      -0.0097     0.0159     0.0080 
# 2019/02/01 00:40:00  58515.027778 2038 434400.000       0.0108       0.0006       0.0140     0.0199     0.0090 
# 2019/02/01 01:00:00  58515.041667 2038 435600.000      -0.0067       0.0116      -0.0096     0.0082     0.0047 
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
        "x",
        "site_vel_x",
        (),
        float,
        "%13.4f",
        13,
        "ECEF [X]",
        "m/s",
        "X-coordinate of station velocity given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "y",
        "site_vel_y",
        (),
        float,
        "%13.4f",
        13,
        "ECEF [Y]",
        "m/s",
        "Y-coordinate of station velocity given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "z",
        "site_vel_z",
        (),
        float,
        "%13.4f",
        13,
        "ECEF [Z]",
        "m/s",
        "Z-coordinate of station velocity given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "sigma_vx",
        "site_vel_sigma_x",
        (),
        float,
        "%11.4f",
        11,
        "SIGMA_X",
        "m/s",
        "Standard deviation of station velocity X-coordinate",
    ),
    WriterField(
        "sigma_vy",
        "site_vel_sigma_y",
        (),
        float,
        "%11.4f",
        11,
        "SIGMA_Y",
        "m/s",
        "Standard deviation of station velocity Y-coordinate",
    ),
    WriterField(
        "sigma_vz",
        "site_vel_sigma_z",
        (),
        float,
        "%11.4f",
        11,
        "SIGMA_Z",
        "m/s",
        "Standard deviation of station velocity Z-coordinate",
    ),
    WriterField(
        "east",
        "site_vel_east",
        (),
        float,
        "%11.4f",
        11,
        "EAST",
        "m/s",
        "East coordinate of station velocity given in topocentric coordinate system",
    ),
    WriterField(
        "north",
        "site_vel_north",
        (),
        float,
        "%11.4f",
        11,
        "NORTH",
        "m/s",
        "North coordinate of station velocity given in topocentric coordinate system",
    ),
    WriterField(
        "up",
        "site_vel_up",
        (),
        float,
        "%11.4f",
        11,
        "UP",
        "m/s",
        "Up coordinate of station velocity given in topocentric coordinate system",
    ),
    WriterField(
        "site_vel_h",
        "site_vel_h",
        (),
        float,
        "%11.4f",
        11,
        "HV",
        "m/s",
        "Horizontal site velocity",
    ),
    WriterField(
        "site_vel_v",
        "site_vel_v",
        (),
        float,
        "%11.4f",
        11,
        "VV",
        "m/s",
        "Vertical site velocity",
    ),
    WriterField(
        "site_vel_3d",
        "site_vel_3d",
        (),
        float,
        "%11.4f",
        11,
        "3D",
        "m/s",
        "3D site velocity",
    ),
    WriterField(
        "c_xx",
        "estimate_cov_site_vel_xx",
        (),
        float,
        "%13.8f",
        13,
        "COV_XX",
        "m**2/s**2",
        "Variance of station velocity X-coordinate",
    ),
    WriterField(
        "c_xy",
        "estimate_cov_site_vel_xy",
        (),
        float,
        "%13.8f",
        13,
        "COV_XY",
        "m**2/s**2",
        "XY-covariance of velocity position",
    ),
    WriterField(
        "c_xz",
        "estimate_cov_site_vel_xz",
        (),
        float,
        "%13.8f",
        13,
        "COV_XZ",
        "m**2/s**2",
        "XZ-covariance of station velocity",
    ),
    WriterField(
        "c_yy",
        "estimate_cov_site_vel_yy",
        (),
        float,
        "%13.8f",
        13,
        "COV_YY",
        "m**2/s**2",
        "Variance of station velocity Y-coordinate",
    ),
    WriterField(
        "c_yz",
        "estimate_cov_site_vel_yz",
        (),
        float,
        "%13.8f",
        13,
        "COV_YZ",
        "m**2/s**2",
        "YZ-covariance of station velocity",
    ),
    WriterField(
        "c_zz",
        "estimate_cov_site_vel_zz",
        (),
        float,
        "%13.8f",
        13,
        "COV_ZZ",
        "m**2/s**2",
        "Variance of station velocity Z-coordinate",
    ),
)


@plugins.register
def gnss_velocity(dset: "Dataset") -> None:
    """Write GNSS velocity results


    Args:
        dset:  A dataset containing the data.
    """
    file_path = config.files.path("output_velocity", file_vars={**dset.vars, **dset.analysis})

    # Add date field to dataset
    if "date" not in dset.fields:
        dset.add_text(
            "date", 
            val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime],
            write_level="detail",
        )

    # Add necessary site velocity fields to dataset    
    if "site_vel_3d" not in dset.fields:
        gnss_velocity_fields(dset)

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

    output_array = np.array(output_list, dtype=[(f.name, f.dtype) for f in FIELDS])

    # Write to disk
    header = get_header(
        FIELDS,
        pgm_version=f"where {where.__version__}",
        run_by=util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else "",
        summary="GNSS velocity results",
    )
    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in FIELDS),
        header=header,
        delimiter="",
        encoding="utf8",
    )
