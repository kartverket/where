"""Write selected RINEX navigation file observations

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
from midgard.math.unit import Unit

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
    WriterField("satellite", "satellite", (), object, "%5s", 5, "SAT", "", "Satellite number"),
    WriterField(
        "iode",
        "iode",
        (),
        float,
        "%6d",
        6,
        "IODE",
        "",
        f"Ephemeris issue of data indicates changes to the broadcast ephemeris:\n"
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
    WriterField("sv_accuracy", "sv_accuracy", (), float, "%8d", 8, "ACC", "", "Satellite vehicle accuracy"),
    WriterField("sv_health", "sv_health", (), float, "%8d", 8, "HLTH", "", "Satellite health status"),
    WriterField("dvs_e1", "dvs_e1", (), float, "%8d", 8, "DVS_E1", "", "Data validity status for E1 signal"),
    WriterField("dvs_e5a", "dvs_e5a", (), float, "%8d", 8, "DVS_E5a", "", "Data validity status for E5a signal"),
    WriterField("dvs_e5b", "dvs_e5b", (), float, "%8d", 8, "DVS_E5b", "", "Data validity status for E5b signal"),
    WriterField("shs_e1", "shs_e1", (), float, "%8d", 8, "SHS_E1", "", "Signal health status for E1 signal"),
    WriterField("shs_e5a", "shs_e5a", (), float, "%8d", 8, "SHS_E5a", "", "Signal health status for E5a signal"),
    WriterField("shs_e5b", "shs_e5b", (), float, "%8d", 8, "SHS_E5b", "", "Signal health status for E5b signal"),
)


@plugins.register
def rinex_nav_writer(dset: "Dataset") -> None:
    """Write selected RINEX navigation file observations

    Args:
        dset:   A dataset containing the data.
    """
    file_path = config.files.path(f"output_rinex_nav", file_vars={**dset.vars, **dset.analysis})

    # Add additional fields used by the writer
    dset.add_text("date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime])
    dset.add_text(
        "time_gpsweek",
        val=[f"{t.gps.gps_ws.week:04.0f}{t.gps.gps_ws.day:1.0f}:{t.gps.gps_ws.seconds:06.0f}" for t in dset.time],
    )
    dset.add_text(
        "trans_time_gpsweek",
        val=[
            f"{t.gps.gps_ws.week:04.0f}{t.gps.gps_ws.day:1.0f}:{t.gps.gps_ws.seconds:06.0f}"
            for t in dset.transmission_time
        ],
    )
    dset.add_text(
        "toe_gpsweek",
        val=[f"{t.gps.gps_ws.week:04.0f}{t.gps.gps_ws.day:1.0f}:{t.gps.gps_ws.seconds:06.0f}" for t in dset.toe],
    )

    dset.add_float("diff_trans_toe", val=(dset.transmission_time.mjd - dset.toe.mjd) * Unit.day2second)

    # Select fields available in Dataset (e.g. DVS and SHS fields are only given for Galileo)
    fields = get_existing_fields(dset, FIELDS)

    # List epochs ordered by satellites
    idx = np.concatenate([np.where(dset.filter(satellite=s))[0] for s in dset.unique("satellite")])

    # Put together fields in an array as specified by the fields-tuple
    output_list = list(zip(*(get_field(dset, f.field, f.attrs, f.unit) for f in fields)))
    output_array = np.array(output_list, dtype=[(f.name, f.dtype) for f in fields])[idx]

    # Write to disk
    header = get_header(
        fields,
        pgm_version=f"where {where.__version__}",
        run_by=util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else "",
        summary="RINEX navigation file analysis results",
    )
    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in fields),
        header=header,
        delimiter="",
        encoding="utf8",
    )
