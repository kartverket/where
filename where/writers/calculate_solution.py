"""Write calculate solution results

Description:
------------


"""
# Standard library imports
from collections import namedtuple

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.writers._writers import get_existing_fields, get_field, get_header

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
# #
# #  SAT               EPOCH    AZI   ELEV          RANGE        SAT_CLK    TROPO  REL_CLK     IONO      TGD            RES
# #      YYYY/MM/DD hh:mm:ss    deg    deg          meter          meter    meter    meter    meter    meter          meter
# # _______________________________________________________________________________________________________________________
#   E02  2019/02/01 00:00:00  -95.2   42.3   24946173.069     -13692.221    3.477    0.180    0.468   -2.303         -1.350
#   E04  2019/02/01 00:00:00  110.7   23.9   26442907.456      63488.395    5.755   -0.037    1.015   -2.583         -0.978
#   E11  2019/02/01 00:00:00   94.2   71.3   23507989.761   -1714889.338    2.473    0.213    0.405   -6.003         -0.520
#   E12  2019/02/01 00:00:00   84.4   17.7   27053455.379   -1819105.131    7.619   -0.066    0.949   -4.886         -0.920
# ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----+----9----+----0----+----1------
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
        "EPOCH",
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
    WriterField(
        "observation", "observation", (), float, "%15.3f", 15, "OBS", "meter", "Observation used in calculation stage"
    ),
    WriterField(
        "range", "delay", ("gnss_range",), float, "%15.3f", 15, "RANGE", "meter", "Station-satellite distance"
    ),
    WriterField(
        "satellite_clock",
        "delay",
        ("gnss_satellite_clock",),
        float,
        "%15.3f",
        15,
        "SAT_CLK",
        "meter",
        "Satellite clock correction",
    ),
    WriterField(
        "troposphere", "delay", ("troposphere_radio",), float, "%9.3f", 9, "TROPO", "meter", "Troposphere delay"
    ),
    WriterField(
        "relativistic_clock",
        "delay",
        ("gnss_relativistic_clock",),
        float,
        "%9.3f",
        9,
        "REL_CLK",
        "meter",
        "Relativistic clock effect due to orbit eccentricity",
    ),
    WriterField("ionosphere", "delay", ("gnss_ionosphere",), float, "%9.3f", 9, "IONO", "meter", "Ionosphere delay"),
    WriterField(
        "satellite_phase_center_offset",
        "delay",
        ("gnss_satellite_phase_center_offset",),
        float,
        "%9.3f",
        9,
        "SAT_OFF",
        "meter",
        "Satellite clock correction",
    ),
    WriterField(
        "total_group_delay",
        "delay",
        ("gnss_total_group_delay",),
        float,
        "%9.3f",
        9,
        "TGD",
        "meter",
        "Total group delay",
    ),
    WriterField(
        "atm_tides_x",
        "site",
        ("atm_tides", "trs", "x"),
        float,
        "%9.3f",
        9,
        "ATMO_X",
        "meter",
        "X-coordinate of site displacement due to atmospheric loading tides related to TRF",
    ),
    WriterField(
        "atm_tides_x",
        "site",
        ("atm_tides", "trs", "y"),
        float,
        "%9.3f",
        9,
        "ATMO_Y",
        "meter",
        "Y-coordinate of site displacement due to atmospheric loading tides related to TRF",
    ),
    WriterField(
        "atm_tides_z",
        "site",
        ("atm_tides", "trs", "z"),
        float,
        "%9.3f",
        9,
        "ATMO_Z",
        "meter",
        "Z-coordinate of site displacement due to atmospheric loading tides related to TRF",
    ),
    WriterField(
        "nt_atm_loading_x",
        "site",
        ("nt_atm_loading", "trs", "x"),
        float,
        "%9.3f",
        9,
        "N_ATMO_X",
        "meter",
        "X-coordinate of site displacement due to non-tidal atmospheric loading related to TRF",
    ),
    WriterField(
        "nt_atm_loading_y",
        "site",
        ("nt_atm_loading", "trs", "y"),
        float,
        "%9.3f",
        9,
        "N_ATMO_Y",
        "meter",
        "Y-coordinate of site displacement due to non-tidal atmospheric loading related to TRF",
    ),
    WriterField(
        "nt_atm_loading_z",
        "site",
        ("nt_atm_loading", "trs", "z"),
        float,
        "%9.3f",
        9,
        "N_ATMO_Z",
        "meter",
        "Z-coordinate of site displacement due to non-tidal atmospheric loading related to TRF",
    ),
    WriterField(
        "ocean_ptides_x",
        "site",
        ("ocean_ptides", "trs", "x"),
        float,
        "%9.3f",
        9,
        "OCN_P_X",
        "meter",
        "X-coordinate of site displacement due to ocean loading tides related to TRF",
    ),
    WriterField(
        "ocean_ptides_y",
        "site",
        ("ocean_ptides", "trs", "y"),
        float,
        "%9.3f",
        9,
        "OCN_P_Y",
        "meter",
        "Y-coordinate of site displacement due to ocean loading tides related to TRF",
    ),
    WriterField(
        "ocean_ptides_z",
        "site",
        ("ocean_ptides", "trs", "z"),
        float,
        "%9.3f",
        9,
        "OCN_P_Z",
        "meter",
        "Z-coordinate of site displacement due to ocean loading tides related to TRF",
    ),
    WriterField(
        "ocean_tides_x",
        "site",
        ("ocean_tides", "trs", "x"),
        float,
        "%9.3f",
        9,
        "OCN_X",
        "meter",
        "X-coordinate of site displacement due to ocean pole tides related to TRF",
    ),
    WriterField(
        "ocean_tides_y",
        "site",
        ("ocean_tides", "trs", "y"),
        float,
        "%9.3f",
        9,
        "OCN_Y",
        "meter",
        "Y-coordinate of site displacement due to ocean pole tides related to TRF",
    ),
    WriterField(
        "ocean_tides_z",
        "site",
        ("ocean_tides", "trs", "z"),
        float,
        "%9.3f",
        9,
        "OCN_Z",
        "meter",
        "Z-coordinate of site displacement due to ocean pole tides related to TRF",
    ),
    WriterField(
        "solid_ptides_x",
        "site",
        ("solid_ptides", "trs", "x"),
        float,
        "%9.3f",
        9,
        "SLD_P_X",
        "meter",
        "X-coordinate fo site displacement due to solid pole tides related to TRF",
    ),
    WriterField(
        "solid_ptides_y",
        "site",
        ("solid_ptides", "trs", "y"),
        float,
        "%9.3f",
        9,
        "SLD_P_Y",
        "meter",
        "Y-coordinate fo site displacement due to solid pole tides related to TRF",
    ),
    WriterField(
        "solid_ptides_z",
        "site",
        ("solid_ptides", "trs", "z"),
        float,
        "%9.3f",
        9,
        "SLD_P_Z",
        "meter",
        "Z-coordinate fo site displacement due to solid pole tides related to TRF",
    ),
    WriterField(
        "solid_tides_x",
        "site",
        ("solid_tides", "trs", "x"),
        float,
        "%9.3f",
        9,
        "SLD_X",
        "meter",
        "X-coordinate of site displacement due to solid tides related to TRF",
    ),
    WriterField(
        "solid_tides_y",
        "site",
        ("solid_tides", "trs", "y"),
        float,
        "%9.3f",
        9,
        "SLD_Y",
        "meter",
        "Y-coordinate of site displacement due to solid tides related to TRF",
    ),
    WriterField(
        "solid_tides_z",
        "site",
        ("solid_tides", "trs", "z"),
        float,
        "%9.3f",
        9,
        "SLD_Z",
        "meter",
        "Z-coordinate of site displacement due to solid tides related to TRF",
    ),
    WriterField("residual_prefit", "residual_prefit", (), float, "%15.3f", 15, "RES", "meter", "Pre-fit residual"),
)


@plugins.register
def calculate_solution(dset: "Dataset") -> None:
    """Write calculate solution results


    Args:
        dset:  A dataset containing the data.
    """
    file_path = config.files.path("output_calculate_solution", file_vars=dset.vars)

    # Add date field to dataset
    if "date" not in dset.fields:
        dset.add_text("date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime])

    # Select fields available in Dataset
    fields = get_existing_fields(dset, FIELDS)

    # Put together fields in an array as specified by the 'dtype' tuple list
    output_list = list(zip(*(get_field(dset, f.field, f.attrs, f.unit) for f in fields)))
    output_array = np.array(output_list, dtype=[(f.name, f.dtype) for f in fields])

    # Write to disk
    header = get_header(
        fields,
        pgm_version=f"where {where.__version__}",
        run_by=util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else "",
        summary="Calculate solutions results",
    )
    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in fields),
        header="\n".join(header),
        delimiter="",
        encoding="utf8",
    )
