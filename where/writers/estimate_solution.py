"""Write estimate solution results

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
from midgard.math.unit import Unit

# Where imports
import where
from where.lib import config
from where.lib import util

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

# Define fields to plot
#
# # PGM: where 0.21.2/midgard 0.3.0  RUN_BY: NMA  DATE: 20190402 194105 UTC
# #
# # MJD          APRIORI                ESTIMATE               SIGMA        PARAM_NAME
# #              meter                  meter                  meter
# # ____________________________________________________________________________________________
#  58150.000000  0.000000000000000e+00 -4.654045691822484e+00  1.747e-01  gnss_rcv_clock
#  58150.000000  3.275756762300001e+06  3.275759445327456e+06  1.485e-01  gnss_site_pos-x
#  58150.000000  3.211111394999998e+05  3.211107445392029e+05  1.038e-01  gnss_site_pos-y
#  58150.000000  5.445046647700000e+06  5.445048372422102e+06  2.453e-01  gnss_site_pos-z
FIELDS = (
    WriterField("date", (), object, "%21s", 19, "DATE", "YYYY/MM/DD hh:mm:ss"),
    WriterField("mjd", (), float, "%13.6f", 13, "MJD", ""),
    WriterField("gpsweek", ("gps", "gpsweek"), int, "%5d", 5, "WEEK", ""),
    WriterField("gpssec", ("gps", "gpssec"), float, "%11.3f", 11, "GPSSEC", "second"),
    WriterField("apriori", (), float, "%23.15e", 23, "APRIORI", "meter"),
    WriterField("estimate", (), float, "%23.15e", 23, "ESTIMATE", "meter"),
    WriterField("sigma", (), float, "%11.3e", 11, "SIGMA", "meter"),
    WriterField("empty", (), str, "%2s", 2, "", ""),  # get spaces between 'sigma' and 'param_name' column
    WriterField("param_name", (), object, "%-20s", 20, "PARAM_NAME", ""),
)


@plugins.register
def estimate_solution(dset: "Dataset") -> None:
    """Write estimate solution results

    Args:
        dset:  A dataset containing the data.
    """
    file_path = config.files.path("output_estimate_solution", file_vars=dset.vars)

    # Add date field to dataset
    if "date" not in dset.fields:
        dset.add_text("date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime])

    # Epochwise estimation or over whole time period
    if config.tech.estimate_epochwise.bool:

        output_array = np.array([])
        for epoch in sorted(set(dset.time.gps.mjd)):
            idx = dset.time.gps.mjd == epoch

            # Append current epoch solution to final output solution for each estimated parameter
            epoch_array = _get_epoch(dset, idx)
            output_array = np.concatenate((output_array, epoch_array), axis=0) if output_array.size else epoch_array
    else:
        # Get solution for first observation
        idx = np.squeeze(np.array(np.nonzero(dset.time.gps.mjd)) == 0)  # first observation -> TODO: Better solution?
        output_array = _get_epoch(dset, idx)

    # Write to disk
    header = [
        _get_header(dset),
        "".join(f"{f.header:>{f.width}s}" for f in FIELDS),
        "".join(f"{f.unit:>{f.width}s}" for f in FIELDS),
        "_" * sum([f.width for f in FIELDS]),
    ]
    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in FIELDS),
        header="\n".join(header),
        delimiter="",
        encoding="utf8",
    )


def _get_epoch(dset: "Dataset", idx: np.ndarray) -> np.ndarray:
    """Get estimation solution for given epoch indices

    Args:
        dset:   A dataset containing the data.
        idx:    Indices for epoch.

    Returns:
        Estimation solution for given epoch indices
    """
    epoch_array = np.array([])

    # Put together fields in an array as specified by the 'dtype' tuple list
    dtype = [(f.field, f.dtype) for f in FIELDS]

    for param_name in dset.meta["param_names"]:

        # Save defined estimation fields epochwise
        values = [
            tuple(
                (
                    dset.date[idx][0],
                    dset.time.gps.mjd[idx][0],
                    dset.time.gps.gps_ws.week[idx][0],
                    dset.time.gps.gps_ws.seconds[idx][0],
                    dset[f"estimate_apriori_{param_name}"][idx][0],
                    dset[f"estimate_{param_name}"][idx][0],
                    dset[f"estimate_sigma_{param_name}"][idx][0],
                    "",
                    param_name,
                )
            )
        ]

        # Append each parameter solution to epoch solution
        param_values = np.array(values, dtype=dtype)
        epoch_array = np.concatenate((epoch_array, param_values), axis=0) if epoch_array.size else param_values

    return epoch_array


def _get_header(dset: "Dataset") -> str:
    """Get header

    Args:
        dset:  A dataset containing the data.

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
DATE        YYYY/MM/DD hh:mm:ss   Date in format year, month, day and hour, minute and second
MJD                               Modified Julian Day
WEEK                              GPS week
GPSSEC      second                GPS seconds
APRIORI     meter                 Apriori value parameter
ESTIMATE    meter                 Estimated parameter value
SIGMA       meter                 Standard deviation of estimated parameter value
PARAM_NAME                        Parameter name:
                                        gnss_rcv_clock   - GNSS receiver clock
                                        gnss_site_pos-x  - X-coordinate of site position
                                        gnss_site_pos-x  - Y-coordinate of site position
                                        gnss_site_pos-x  - Z-coordinate of site position

"""
    )
    return header
