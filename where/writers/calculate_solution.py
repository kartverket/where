"""Write calculate solution results

Description:
------------


"""
# Standard library imports
from collections import namedtuple
from datetime import datetime
from typing import List, Tuple

# External library imports
import numpy as np

# Midgard imports
import midgard
from midgard.dev import plugins

# Where imports
import where
from where.lib import config
from where.lib.unit import Unit
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
# #  SAT               EPOCH           MJD WEEK     GPSSEC    AZI   ELEV            OBS          RANGE        SAT_CLK
# #      YYYY/MM/DD hh:mm:ss                        second    deg    deg          meter          meter          meter
# # __________________________________________________________________________________________________________________
#   E13  2019/07/01 00:00:00  58665.000000 2060  86400.000  241.1    5.4   28190341.078   28307032.202    -116713.398
#   E30  2019/07/01 00:00:00  58665.000000 2060  86400.000  184.3    5.7   26839718.930   28262933.516   -1423236.207
#   E25  2019/07/01 00:00:00  58665.000000 2060  86400.000   67.5   44.6   24267086.023   24793731.055    -526646.818
# ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----+----9----+----0----+----1--
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
    # TODO: Site model corrections are given for XYZ. This has to be handled.
    # WriterField("atmospheric_tides", "site", ("atmospheric_tides",), float, "%9.3f", 9, "ATMO_TDS", "meter", "Site displacement due to atmospheric loading tides"),
    # WriterField("non_tidal_atmospheric_loading", "site", ("non_tidal_atmospheric_loading",), float, "%9.3f", 9, "N_ATMO_TDS", "meter", "Site displacement due to non-tidal atmospheric loading"),
    # WriterField("ocean_pole_tides", "site", ("ocean_pole_tides",), float, "%9.3f", 9, "OCN_POL_TDS", "meter", "Site displacement due to ocean loading tides."),
    # WriterField("ocean_tides", "site", ("ocean_tides",), float, "%9.3f", 9, "OCN_TDS", "meter", "Site displacement due to ocean pole tides."),
    # WriterField("solid_pole_tides", "site", ("solid_pole_tides",), float, "%9.3f", 9, "SLD_POL_TDS", "meter", "Site displacement due to solid pole tides."),
    # WriterField("solid_tides", "site", ("solid_tides",), float, "%9.3f", 9, "SLD_TDS", "meter", "Site displacement due to solid tides."),
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
    fields = []
    for f in FIELDS:
        field_attrs = f.field if len(f.attrs) == 0 else f"{f.field}.{'.'.join(f.attrs)}"
        if field_attrs in dset.fields:
            fields.append(f)

    # Put together fields in an array as specified by the 'dtype' tuple list
    output_list = list(zip(*(_get_field(dset, f.field, f.attrs, f.unit) for f in fields)))
    output_array = np.array(output_list, dtype=[(f.name, f.dtype) for f in fields])

    # Write to disk
    header = [
        _get_header(fields),
        "".join(f"{f.header:>{f.width}s}" for f in fields),
        "".join(f"{f.unit:>{f.width}s}" for f in fields),
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


def _get_field(dset: "Dataset", field: "str", attrs: Tuple[str], unit: "str") -> np.ndarray:
    """Get field values of a Dataset specified by the field attributes

    If necessary the unit of the data fields are corrected to the defined 'output' unit.

    Args:
        dset:     Dataset, a dataset containing the data.
        field:    Field name.
        attrs:    Field attributes (e.g. for Time object: (<scale>, <time format>)).
        unit:     Unit used for output.

    Returns:
        Array with Dataset field values
    """
    f = dset[field]
    for attr in attrs:
        f = getattr(f, attr)

    # Determine output 'unit'
    # +TODO: Does not work for all fields, because dset.unit() does not except 'time.gps.mjd'.
    if unit.startswith("deg"):
        field_attrs = field if len(attrs) == 0 else f"{field}.{'.'.join(attrs)}"
        f = f * getattr(Unit, f"{dset.unit(field_attrs)[0]}2{unit}")
    # -TODO

    return f


def _get_header(fields: List["str"]) -> str:
    """Get header

    Args:
        fields:  List with fields to write.

    Returns:
        Header lines
    """
    pgm = "where " + where.__version__ + "/midgard " + midgard.__version__
    run_by = util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else ""
    file_created = datetime.utcnow().strftime("%Y%m%d %H%M%S") + " UTC"
    header = "PGM: {:s}  RUN_BY: {:s}  DATE: {:s}\n".format(pgm, run_by, file_created)
    header = (
        header
        + """

HEADER      UNIT                  DESCRIPTION
______________________________________________________________________________________________________________________
"""
    )
    for f in fields:
        header = f"{header}{f.header:11s} {f.unit:22s} {f.description}\n"

    return header
