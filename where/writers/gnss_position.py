"""Write GNSS position results

Description:
------------


"""
# Standard library imports
from collections import namedtuple

# External library imports
import numpy as np

# Midgard imports
from midgard.data import position
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
# # PGM: where 0.21.2/midgard 0.3.0  RUN_BY: NMA  DATE: 20190604 135301 UTC
# #
# #               EPOCH           MJD WEEK     GPSSEC     ECEF [X]     ECEF [Y]     ECEF [Z]     LATITUDE
# # YYYY/MM/DD hh:mm:ss                        second            m            m            m          deg
# # ___________________________________________________________________________________________________________
#   2018/02/01 00:00:00  58150.000000 1986 345600.000 3275756.9411  321112.5306 5445048.1920  59.01771380   ...
#   2018/02/01 00:05:00  58150.003472 1986 345900.000 3275756.9004  321112.5296 5445048.4119  59.01771513   ...
#   2018/02/01 00:10:00  58150.006944 1986 346200.000 3275757.1458  321112.6499 5445047.0430  59.01770683   ...
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
        "site_pos",
        ("trs", "x"),
        float,
        "%13.4f",
        13,
        "ECEF [X]",
        "meter",
        "X-coordinate of station position given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "y",
        "site_pos",
        ("trs", "y"),
        float,
        "%13.4f",
        13,
        "ECEF [Y]",
        "meter",
        "Y-coordinate of station position given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "z",
        "site_pos",
        ("trs", "z"),
        float,
        "%13.4f",
        13,
        "ECEF [Z]",
        "meter",
        "Z-coordinate of station position given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "sigma_x",
        "site_pos_sigma_x",
        (),
        float,
        "%11.4f",
        11,
        "SIGMA_X",
        "meter",
        "Standard deviation of station position X-coordinate",
    ),
    WriterField(
        "sigma_y",
        "site_pos_sigma_y",
        (),
        float,
        "%11.4f",
        11,
        "SIGMA_Y",
        "meter",
        "Standard deviation of station position Y-coordinate",
    ),
    WriterField(
        "sigma_z",
        "site_pos_sigma_z",
        (),
        float,
        "%11.4f",
        11,
        "SIGMA_Z",
        "meter",
        "Standard deviation of station position Z-coordinate",
    ),
    WriterField(
        "lat",
        "site_pos",
        ("llh", "lat"),
        float,
        "%13.8f",
        13,
        "LATITUDE",
        "degree",
        "Latitude coordinate of station position given in ellipsiodal reference frame",
    ),
    WriterField(
        "lon",
        "site_pos",
        ("llh", "lon"),
        float,
        "%13.8f",
        13,
        "LONGITUDE",
        "degree",
        "Longitude coordinate of station position given in ellipsiodal reference frame",
    ),
    WriterField(
        "h",
        "site_pos",
        ("llh", "height"),
        float,
        "%11.4f",
        11,
        "HEIGHT",
        "meter",
        "Height coordinate of station position given in ellipsiodal reference frame",
    ),
    WriterField(
        "east",
        "site_pos_vs_ref_east",
        (),
        float,
        "%11.4f",
        11,
        "EAST",
        "meter",
        "East coordinate of difference between station position and reference position (e.g."
        "apriori station coordinate) given in topocentric coordinate system",
    ),
    WriterField(
        "north",
        "site_pos_vs_ref_north",
        (),
        float,
        "%11.4f",
        11,
        "NORTH",
        "meter",
        "North coordinate of difference between station position and reference position (e.g."
        "apriori station coordinate) given in topocentric coordinate system",
    ),
    WriterField(
        "up",
        "site_pos_vs_ref_up",
        (),
        float,
        "%11.4f",
        11,
        "UP",
        "meter",
        "Up coordinate of difference between station position and reference position (e.g."
        "apriori station coordinate) given in topocentric coordinate system",
    ),
    WriterField(
        "hpe",
        "hpe",
        (),
        float,
        "%11.4f",
        11,
        "HPE",
        "meter",
        "Horizontal Position Error of station position vs. reference position",
    ),
    WriterField(
        "vpe",
        "vpe",
        (),
        float,
        "%11.4f",
        11,
        "VPE",
        "meter",
        "Vertical Position Error of station position vs. reference position",
    ),
    WriterField(
        "c_xx",
        "estimate_cov_site_pos_xx",
        (),
        float,
        "%13.8f",
        13,
        "COV_XX",
        "meter**2",
        "Variance of station position X-coordinate",
    ),
    WriterField(
        "c_xy",
        "estimate_cov_site_pos_xy",
        (),
        float,
        "%13.8f",
        13,
        "COV_XY",
        "meter**2",
        "XY-covariance of station position",
    ),
    WriterField(
        "c_xz",
        "estimate_cov_site_pos_xz",
        (),
        float,
        "%13.8f",
        13,
        "COV_XZ",
        "meter**2",
        "XZ-covariance of station position",
    ),
    WriterField(
        "c_yy",
        "estimate_cov_site_pos_yy",
        (),
        float,
        "%13.8f",
        13,
        "COV_YY",
        "meter**2",
        "Variance of station position Y-coordinate",
    ),
    WriterField(
        "c_yz",
        "estimate_cov_site_pos_yz",
        (),
        float,
        "%13.8f",
        13,
        "COV_YZ",
        "meter**2",
        "YZ-covariance of station position",
    ),
    WriterField(
        "c_zz",
        "estimate_cov_site_pos_zz",
        (),
        float,
        "%13.8f",
        13,
        "COV_ZZ",
        "meter**2",
        "Variance of station position Z-coordinate",
    ),
)


@plugins.register
def gnss_position(dset: "Dataset") -> None:
    """Write GNSS position results


    Args:
        dset:  A dataset containing the data.
    """
    file_path = config.files.path("output_position", file_vars={**dset.vars, **dset.analysis})

    # Add date field to dataset
    if "date" not in dset.fields:
        dset.add_text("date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime], write_level="detail")

    # Add ENU position to dataset
    ref_pos = position.Position(
        val=np.array([dset.meta["pos_x"], dset.meta["pos_y"], dset.meta["pos_z"]]), system="trs"
    )
    enu = (dset.site_pos.trs.pos - ref_pos).enu
    dset.add_float("site_pos_vs_ref_east", val=enu.east, unit="meter", write_level="detail")
    dset.add_float("site_pos_vs_ref_north", val=enu.north, unit="meter", write_level="detail")
    dset.add_float("site_pos_vs_ref_up", val=enu.up, unit="meter", write_level="detail")

    # Add HPE and VPE to dataset
    dset.add_float("hpe", val=np.sqrt(enu.east ** 2 + enu.north ** 2), unit="meter", write_level="operational")
    dset.add_float("vpe", val=np.absolute(enu.up), unit="meter", write_level="operational")

    # Add standard deviation of site position coordinates
    dset.add_float(
            "site_pos_sigma_x", 
            val=np.sqrt(dset.estimate_cov_site_pos_xx), 
            unit="meter",
            write_level="detail",
    )
    dset.add_float(
            "site_pos_sigma_y", 
            val=np.sqrt(dset.estimate_cov_site_pos_yy), 
            unit="meter",
            write_level="detail",
    )
    dset.add_float(
            "site_pos_sigma_z", 
            val=np.sqrt(dset.estimate_cov_site_pos_zz), 
            unit="meter",
            write_level="detail",
    )

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
        summary="GNSS position results",
    )
    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in FIELDS),
        header=header,
        delimiter="",
        encoding="utf8",
    )
