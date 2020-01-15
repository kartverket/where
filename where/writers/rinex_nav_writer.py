"""Write selected RINEX navigation file observations

Description:
------------


"""
# Standard library imports
from collections import namedtuple
from datetime import datetime
from typing import Tuple

# External library imports
import numpy as np

# Midgard imports
import midgard
from midgard.dev import console
from midgard.dev import plugins

# Where imports
import where
from where.lib import config
from where.lib.unit import Unit
from where.lib import util
from where import pipelines

WriterField = namedtuple("WriterField", ["field", "attrs", "dtype", "format", "width", "header", "unit"])
WriterField.__new__.__defaults__ = (None,) * len(WriterField._fields)
WriterField.__doc__ = """A convenience class for defining a output field for the writer

    Args:
        field (str):             Dataset field name
        attrs (Tuple[str]):      Field attributes
        dtype (Numpy dtype):     Type of field
        format (str):            Format string
        width (int):             Width of header information
        header (str):            Header information
        unit (str):              Unit of field
    """


@plugins.register
def rinex_nav_writer(dset: "Dataset") -> None:
    """Write selected RINEX navigation file observations

    Args:
        dset:   A dataset containing the data.
    """
    write_level = config.tech.get("write_level", default="operational").as_enum("write_level")

    fields = (
        WriterField("time_date", (), object, "%21s", 19, "EPOCH", "YYYY/MM/DD hh:mm:ss"),
        WriterField("time", ("gps", "mjd"), float, "%14.6f", 14, "", "mjd"),
        WriterField("time_gpsweek", (), object, "%15s", 15, "", "wwwwd:ssssss"),
        WriterField("satellite", (), object, "%5s", 5, "SAT", ""),
        WriterField("iode", (), float, "%6d", 6, "IODE", ""),
        WriterField("trans_time_gpsweek", (), object, "%15s", 15, "TRANS_TIME", "wwwwd:ssssss"),
        WriterField("toe_gpsweek", (), object, "%15s", 15, "TOE", "wwwwd:ssssss"),
        WriterField("diff_trans_toe", (), float, "%8d", 8, "TM-TOE", "second"),
        WriterField("sv_accuracy", (), float, "%8d", 8, "ACC", ""),
        WriterField("sv_health", (), float, "%8d", 8, "HLTH", ""),
    )

    additional_galileo_fields = {
        "dvs_e1": WriterField("dvs_e1", (), float, "%8d", 8, "DVS_E1", ""),
        "dvs_e5a": WriterField("dvs_e5a", (), float, "%8d", 8, "DVS_E5a", ""),
        "dvs_e5b": WriterField("dvs_e5b", (), float, "%8d", 8, "DVS_E5b", ""),
        "shs_e1": WriterField("shs_e1", (), float, "%8d", 8, "SHS_E1", ""),
        "shs_e5a": WriterField("shs_e5a", (), float, "%8d", 8, "SHS_E5a", ""),
        "shs_e5b": WriterField("shs_e5b", (), float, "%8d", 8, "SHS_E5b", ""),
    }

    # Add additional fields used by the writer
    dset.add_text("time_date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime])
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

    # Add additional Galileo fields if available
    for field in sorted(additional_galileo_fields.keys()):
        if field in dset.fields:
            fields = list(fields)
            fields.append(additional_galileo_fields[field])
            fields = tuple(fields)
            # fields += additional_galileo_fields[field]

    # List epochs ordered by satellites
    idx = np.concatenate([np.where(dset.filter(satellite=s))[0] for s in dset.unique("satellite")])

    # Put together fields in an array as specified by the fields-tuple
    output_list = list(zip(*(_get_field(dset, f.field, f.attrs) for f in fields)))
    output_array = np.array(output_list, dtype=[(f.field, f.dtype) for f in fields])[idx]

    # Write to disk
    file_path = config.files.path(f"output_rinex_nav", file_vars=dset.vars)
    header = [
        _get_header(dset),
        "".join(f"{f.header:>{f.width}s}" for f in fields),
        "".join(f"{f.unit if f.unit is not None else dset.unit(f.field):>{f.width}s}" for f in fields),
        "_" * sum([f.width for f in fields]),
    ]
    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in fields),
        header="\n".join(header),
        delimiter="",
        encoding="utf8",
    )


def _get_field(dset: "Dataset", field: "str", attrs: Tuple[str]) -> np.ndarray:
    """Get field values of a Dataset specified by the field attributes

    Args:
        dset:     Dataset, a dataset containing the data.
        field:    Field name.
        attrs:    Field attributes (e.g. for Time object: (<scale>, <time format>)).

    Returns:
        Array with Dataset field values
    """
    f = dset[field]
    for attr in attrs:
        f = getattr(f, attr)
    return f


def _get_header(dset: "Dataset") -> str:
    """Get header

    Args:
        dset:   A dataset containing the data.

    Returns:
        Header lines
    """

    pgm = "where " + where.__version__ + "/midgard " + midgard.__version__
    run_by = util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else ""
    file_created = datetime.utcnow().strftime("%Y%m%d %H%M%S") + " UTC"
    header = "PGM: {:s}  RUN_BY: {:s}  DATE: {:s}\n\n".format(pgm, run_by, file_created)
    header += "RINEX NAVIGATION ANALYSIS CONFIGURATION\n\n"

    # RINEX navigation configuration
    header += str(config.tech.as_str(key_width=25, width=70, only_used=True)) + "\n\n\n"
    header += _get_paths()

    header = (
        header
        + """

HEADER      UNIT                  DESCRIPTION
______________________________________________________________________________________________________________________
DATE        YYYY/MM/DD hh:mm:ss   Date in format year, month, day and hour, minute and second
MJD                               Modified Julian Day
WEEK                              GPS week and day
GPSSEC      second                GPS seconds
SAT                               Satellite number
IODE                              Ephemeris issue of data indicates changes to the broadcast ephemeris:
                                       - GPS:     Ephemeris issue of data (IODE), which is set equal to IODC
                                       - Galileo: Issue of Data of the NAV batch (IODnav)
                                       - QZSS:    Ephemeris issue of data (IODE)
                                       - BeiDou:  Age of Data Ephemeris (AODE)
                                       - IRNSS:   Issue of Data, Ephemeris and Clock (IODEC)
TRANS_TIME  WWWWD:SSSSS           Transmission time (receiver reception time)
TOE         WWWWD:SSSSS           Time of ephemeris
TM-TOE      second                Difference between transmission time and time of ephemeris
ACC                               Satellite vehicle accuracy
HLTH                              Satellite health status
"""
    )
    additional_galileo_header_info = {
        "dvs_e1": "DVS_E1                            Data validity status for E1 signal\n",
        "dvs_e5a": "DVS_E5a                           Data validity status for E5a signal\n",
        "dvs_e5b": "DVS_E5b                           Data validity status for E5b signal\n",
        "shs_e1": "SHS_E1                            Signal health status for E1 signal\n",
        "shs_e5a": "SHS_E5a                           Signal health status for E5a signal\n",
        "shs_e5b": "SHS_E5b                           Signal health status for E5b signal\n",
    }

    # Add additional Galileo fields if available
    for field in sorted(additional_galileo_header_info.keys()):
        if field in dset.fields:
            header += additional_galileo_header_info[field]

    return header


def _get_paths() -> str:
    """Get file paths of used files

    Returns:
        Header file path lines
    """
    lines = "[file_paths]\n"
    key_width = 25
    fill_args = dict(width=120, hanging=key_width + 3, break_long_words=False, break_on_hyphens=False)
    path_def = {"Broadcast orbit": "gnss_rinex_nav_."}

    for name, file_key in sorted(path_def.items()):
        path = ", ".join(str(p) for p in sorted(pipelines.paths(file_key)))
        lines += console.fill(f"{name:<{key_width}} = {path}", **fill_args)
        lines += "\n"

    lines += "\n\n"
    return lines
