"""Write GNSS satellite position results

Description:
------------


"""
# Standard library imports
from collections import namedtuple

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit
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
# # PGM: where 1.0.4/midgard 1.1.1  RUN_BY: NMA  DATE: 20201001 143801 UTC
# # DESCRIPTION: GNSS satellite position results
# #    
# #                DATE           MJD WEEK     GPSSEC  SAT        ECEF [X]        ECEF [Y]        ECEF [Z]    ...
# # YYYY/MM/DD hh:mm:ss                        second                meter           meter           meter    ... 
# # ______________________________________________________________________________________________________________
#    
#   2019/02/01 00:00:00  58515.000000 2038 432000.000  E02   15095082.6158  -16985925.1555   18975783.7802    ...
#   2019/02/01 00:00:00  58515.000000 2038 432000.000  E04   13831647.1962   24089264.5694   10227973.9704    ...    
#   2019/02/01 00:00:00  58515.000000 2038 432000.000  E11   14423097.6220    8967728.5570   24249267.0732    ...  
#   2019/02/01 00:00:00  58515.000000 2038 432000.000  E12    2848971.0942   26052257.7493   13788700.2981    ...    
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
    WriterField("satellite", "satellite", (), object, "%5s", 5, "SAT", " ", "Satellite number"),
    WriterField(
        "x",
        "sat_posvel",
        ("trs", "x"),
        float,
        "%16.4f",
        16,
        "ECEF [X]",
        "meter",
        "X-coordinate of satellite position given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "y",
        "sat_posvel",
        ("trs", "y"),
        float,
        "%16.4f",
        16,
        "ECEF [Y]",
        "meter",
        "Y-coordinate of satellite position given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "z",
        "sat_posvel",
        ("trs", "z"),
        float,
        "%16.4f",
        16,
        "ECEF [Z]",
        "meter",
        "Z-coordinate of satellite position given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "vx",
        "sat_posvel",
        ("trs", "vx"),
        float,
        "%13.4f",
        13,
        "ECEF [VX]",
        "meter/second",
        "X-coordinate of satellite velocity given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "vy",
        "sat_posvel",
        ("trs", "vy"),
        float,
        "%13.4f",
        13,
        "ECEF [VY]",
        "meter/second",
        "Y-coordinate of satellite velocity given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "vz",
        "sat_posvel",
        ("trs", "vz"),
        float,
        "%13.4f",
        13,
        "ECEF [VZ]",
        "meter/second",
        "Z-coordinate of satellite velocity given in Earth-Centered Earth-Fixed cartesian coordinate system",
    ),
    WriterField(
        "used_iode",
        "used_iode",
        (),
        float,
        "%6d",
        6,
        "IODE",
        "",
        f"Used ephemeris issue of data indicates changes to the broadcast ephemeris:\n"
        f"""
{'': >38}- GPS:     Ephemeris issue of data (IODE), which is set equal to IODC
{'': >38}- Galileo: Issue of Data of the NAV batch (IODnav)
{'': >38}- QZSS:    Ephemeris issue of data (IODE)
{'': >38}- BeiDou:  Age of Data Ephemeris (AODE)
""",
    ),
    WriterField(
        "trans_time_gpsweek",
        "trans_time_gpsweek",
        (),
        object,
        "%15s",
        15,
        "TRANS_TIME",
        "wwwwd:ssssss",
        "Transmission time (receiver reception time)",
    ),
    WriterField("toe_gpsweek", "toe_gpsweek", (), object, "%15s", 15, "TOE", "wwwwd:ssssss", "Time of ephemeris"),
    WriterField(
        "diff_trans_toe",
        "diff_trans_toe",
        (),
        float,
        "%8d",
        8,
        "TM-TOE",
        "second",
        "Difference between transmission time and time of ephemeris",
    ),
    WriterField(
        "age_of_ephemeris",
        "age_of_ephemeris",
        (),
        float,
        "%8d",
        8,
        "T-TOE",
        "second",
        "Age of ephemeris, which is the difference between the observation time and the time of ephemeris " "(ToE)",
    ),
)


@plugins.register
def gnss_satellite_position(dset: "Dataset") -> None:
    """Write GNSS satellite position results


    Args:
        dset:  A dataset containing the data.
    """
    file_path = config.files.path("output_satellite_position", file_vars={**dset.vars, **dset.analysis})

    # Add date field to dataset
    if "date" not in dset.fields:
        dset.add_text("date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime], write_level="detail")
        
    # Add fields in case of broadcast ephemeris
    if config.tech.apriori_orbit.str == "broadcast":
        dset.add_text(
            "trans_time_gpsweek",
            val=[
                f"{t.gps_ws.week:04.0f}{t.gps_ws.day:1.0f}:{t.gps_ws.seconds:06.0f}" for t in dset.used_transmission_time
            ],
            write_level="detail",
        )
        dset.add_text(
            "toe_gpsweek",
            val=[f"{t.gps_ws.week:04.0f}{t.gps_ws.day:1.0f}:{t.gps_ws.seconds:06.0f}" for t in dset.used_toe],
            write_level="detail",
        )
        dset.add_float(
            "diff_trans_toe",
            val=(dset.used_transmission_time.gps.mjd - dset.used_toe.gps.mjd) * Unit.day2second,
            unit="second", 
            write_level="detail",
        )
        dset.add_float(
            "age_of_ephemeris",
            val=(dset.time.gps.mjd - dset.used_toe.gps.mjd) * Unit.day2second,
            unit="second", 
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
        summary="GNSS satellite position results",
    )
    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in fields),
        header=header,
        delimiter="",
        encoding="utf8",
    )
