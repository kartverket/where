"""Write Galileo HAS correction

Description:
------------


"""
# Standard library imports
from collections import namedtuple

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import log
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
# #
# # PGM: where 1.3.0/midgard 1.1.5  RUN_BY: NMA  DATE: 20220524 124125 UTC
# # DESCRIPTION: Galileo HAS correction
# # 
# # 
# # 
# # HEADER         UNIT                   DESCRIPTION
# # _____________________________________________________________________________________________________________________
# # DATE           YYYY/MM/DD hh:mm:ss    Date in format year, month, day and hour, minute and second
# # MJD                                   Modified Julian Day
# # WEEK                                  GPS week
# # GPSSEC         second                 GPS seconds
# # SAT                                   Satellite number
# # BRC_ORB_X      meter                  X-coordinate of broadcast orbits given in Earth-Centered Earth-Fixed cartesian coordinate system
# # BRC_ORB_Y      meter                  Y-coordinate of broadcast orbits given in Earth-Centered Earth-Fixed cartesian coordinate system
# # BRC_ORB_Z      meter                  Z-coordinate of broadcast orbits given in Earth-Centered Earth-Fixed cartesian coordinate system
# # HAS_ORB [X]    meter                  X-coordinate of Galileo HAS orbit correction given in Earth-Centered Earth-Fixed cartesian coordinate system
# # HAS_ORB [Y]    meter                  Y-coordinate of Galileo HAS orbit correction given in Earth-Centered Earth-Fixed cartesian coordinate system
# # HAS_ORB [Z]    meter                  Z-coordinate of Galileo HAS orbit correction given in Earth-Centered Earth-Fixed cartesian coordinate system
# # BRC_CLK        meter                  Broadcast satellite clock correction
# # HAS_CLK        meter                  Galileo HAS clock correction
# # CBIAS_C1X      meter                  Galileo HAS code bias correction for E1-C signal used for C1X observation
# # CBIAS_C7X      meter                  Galileo HAS code bias correction for E5b-Q signal used for C7X observation
# # BRC_IOD                               Issue of Data (IOD) of broadcast navigation messages
# # HAS_IOD                               Reference Issue of Data (IOD), which corresponds to IOD of broadcast navigation messages
# # 
# #                DATE           MJD WEEK     GPSSEC SAT       BRC_ORB_X       BRC_ORB_Y       BRC_ORB_Z  ...
# # YYYY/MM/DD hh:mm:ss                        second               meter           meter           meter  ...
# # __________________________________________________________________________________________________________
# # 
#   2022/03/21 11:55:00  59659.496528 2202 129300.000 E21    8525157.7714  -24230508.3018   14695938.3363  ...
#   2022/03/21 11:55:00  59659.496528 2202 129300.000 E27   16228222.1517   -4919748.2483   24257198.2708  ...
#   2022/03/21 11:55:00  59659.496528 2202 129300.000 E15   25031045.8205   15204514.4336    4256044.5803  ...
#   2022/03/21 11:55:00  59659.496528 2202 129300.000 E30   13895318.4396   18168221.9722   18785574.3303  ...
#   2022/03/21 11:55:00  59659.496528 2202 129300.000 E04  -17903452.0655   -7137106.1475   22474453.1836  ...
#   2022/03/21 12:00:00  59659.500000 2202 129600.000 E21    8535817.8700  -23776543.4725   15413604.4998  ...
#   2022/03/21 12:00:00  59659.500000 2202 129600.000 E27   16278067.6142   -4189489.0259   24361045.7122  ...
# # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----
# #
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
        "Date in format year, month, day and hour, minute and second of observation epoch",
    ),
    WriterField(
        "mjd", 
        "time", 
        ("gps", "mjd"), 
        float, 
        "%14.6f", 
        14, 
        "MJD", 
        "", 
        "Modified Julian Day of observation epoch",
    ),
    WriterField(
        "gpsweek", 
        "time", 
        ("gps", "gps_ws", "week"), 
        int, 
        "%5d", 
        5,
        "WEEK", 
        "", 
        "GPS week of observation epoch",
    ),
    WriterField(
        "gpssec", 
        "time", 
        ("gps", "gps_ws", "seconds"), 
        float, 
        "%11.0f", 
        11, 
        "GPSSEC", 
        "second", 
        "GPS seconds of observation epoch",
    ),
    WriterField(
        "tom_gpssec_orb", 
        "has_time_of_message_orb", 
        ("gps", "gps_ws", "seconds"), 
        float, 
        "%11.0f", 
        11, 
        "TOM_ORB", 
        "second", 
        "GPS seconds from reference time of HAS orbit correction message",
    ),
    WriterField(
        "rec_gpssec_orb", 
        "has_reception_time_of_message_orb", 
        ("gps", "gps_ws", "seconds"), 
        float, 
        "%11.0f", 
        11, 
        "REC_ORB", 
        "second", 
        "GPS seconds from receiver reception time of HAS orbit correction message",
    ),
    WriterField(
        "tom_gpssec_clk", 
        "has_time_of_message_clk", 
        ("gps", "gps_ws", "seconds"), 
        float, 
        "%11.0f", 
        11, 
        "TOM_CLK", 
        "second", 
        "GPS seconds from reference time of HAS clock correction message",
    ),
    WriterField(
        "rec_gpssec_clk", 
        "has_reception_time_of_message_clk", 
        ("gps", "gps_ws", "seconds"), 
        float, 
        "%11.0f", 
        11, 
        "REC_CLK", 
        "second", 
        "GPS seconds from receiver reception time of HAS clock correction message",
    ),
    WriterField(
        "used_iode",
        "used_iode",
        (),
        float,
        "%8.0f",
        8,
        "IOD_NAV",
        "",
        "Issue of Data (IOD) of broadcast navigation messages",
    ),
    WriterField(
        "has_gnssiod_orb",
        "has_gnssiod_orb",
        (),
        float,
        "%8.0f",
        8,
        "IOD_ORB",
        "",
        "Reference Issue of Data (IOD) of HAS orbit correction message, which corresponds to IOD of broadcast navigation messages",
    ),
    WriterField(
        "has_gnssiod_clk",
        "has_gnssiod_clk",
        (),
        float,
        "%8.0f",
        8,
        "IOD_CLK",
        "",
        "Reference Issue of Data (IOD) of HAS clock correction message, which corresponds to IOD of broadcast navigation messages",
    ),
    WriterField(
        "satellite",
        "satellite",
        (),
        object,
        "%5s",
        5,
        "SAT",
        "",
        "Satellite number",
    ),
    WriterField(
        "used_satpos_x",
        "used_satpos_x",
        (),
        float,
        "%16.4f",
        16,
        "BRC_ORB_X",
        "meter",
        "X-coordinate of broadcast orbits given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "used_satpos_y",
        "used_satpos_y",
        (),
        float,
        "%16.4f",
        16,
        "BRC_ORB_Y",
        "meter",
        "Y-coordinate of broadcast orbits given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "used_satpos_z",
        "used_satpos_z",
        (),
        float,
        "%16.4f",
        16,
        "BRC_ORB_Z",
        "meter",
        "Z-coordinate of broadcast orbits given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "x",
        "has_orbit_correction",
        ("trs", "x"),
        float,
        "%13.4f",
        13,
        "HAS_ORB [X]",
        "meter",
        "X-coordinate of Galileo HAS orbit correction given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "y",
        "has_orbit_correction",
        ("trs", "y"),
        float,
        "%13.4f",
        13,
        "HAS_ORB [Y]",
        "meter",
        "Y-coordinate of Galileo HAS orbit correction given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "z",
        "has_orbit_correction",
        ("trs", "z"),
        float,
        "%13.4f",
        13,
        "HAS_ORB [Z]",
        "meter",
        "Z-coordinate of Galileo HAS orbit correction given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "used_gnss_satellite_clock",
        "used_gnss_satellite_clock",
        (),
        float,
        "%14.4f",
        14,
        "BRC_CLK",
        "meter",
        "Broadcast satellite clock correction",
    ),
    WriterField(
        "has_clock_correction",
        "has_clock_correction",
        (),
        float,
        "%13.4f",
        13,
        "HAS_CLK",
        "meter",
        "Galileo HAS clock correction",
    ),
    WriterField(
        "has_code_bias_C1X",
        "has_code_bias",
        ("C1X",),
        float,
        "%13.4f",
        13,
        "CBIAS_C1X",
        "meter",
        "Galileo HAS code bias correction for E1-C signal used for C1X observation",
    ),
    WriterField(
        "has_code_bias_C5X",
        "has_code_bias",
        ("C5X",),
        float,
        "%13.4f",
        13,
        "CBIAS_C5X",
        "meter",
        "Galileo HAS code bias correction for E5a-Q signal used for C5X observation",
    ),
    WriterField(
        "has_code_bias_C7X",
        "has_code_bias",
        ("C7X",),
        float,
        "%13.4f",
        13,
        "CBIAS_C7X",
        "meter",
        "Galileo HAS code bias correction for E5b-Q signal used for C7X observation",
    ),
)


@plugins.register
def gnss_has_correction(dset: "Dataset") -> None:
    """Write Galileo HAS correction


    Args:
        dset:  A dataset containing the data.
    """  
    file_path = config.files.path("output_has_correction", file_vars={**dset.vars, **dset.analysis})
    
    # Skip using this writer, if 
    if not "has_orbit_correction" in dset.fields:
        log.warn(f"No Galileo HAS correction in dataset. Therefore file '{file_path}' is not written.")
        return

    # Add date field to dataset
    if "date" not in dset.fields:
        dset.add_text("date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime], write_level="detail")
    
    # Add original orbit and clock correction
    dset.add_float(
        "used_satpos_x", val=dset.sat_posvel.trs.x - dset.has_orbit_correction.pos.trs.x, 
        unit="meter", 
        write_level="detail",
    )
    dset.add_float(
        "used_satpos_y", val=dset.sat_posvel.trs.y - dset.has_orbit_correction.pos.trs.y, 
        unit="meter", 
        write_level="detail",
    )
    dset.add_float(
        "used_satpos_z", 
        val=dset.sat_posvel.trs.z - dset.has_orbit_correction.pos.trs.z, 
        unit="meter",
        write_level="detail",
    )
    dset.add_float(
        "used_gnss_satellite_clock", 
        val=-(dset.delay.gnss_satellite_clock + dset.has_clock_correction), 
        unit="meter",
        write_level="detail",
    )

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
        summary="Galileo HAS correction",
    )

    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in fields),
        header=header,
        delimiter="",
        encoding="utf8",
    )
