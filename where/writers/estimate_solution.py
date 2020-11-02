"""Write estimate solution results

Description:
------------


"""
# Standard library imports
from collections import namedtuple
from typing import List, Tuple

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.writers._writers import get_header

# Where imports
import where
from where.lib import config
from where.lib import log
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
    WriterField("mjd", "mjd", (), float, "%13.6f", 13, "MJD", "", "Modified Julian Day"),
    WriterField("gpsweek", "gpsweek", ("gps", "gps_ws", "week"), int, "%5d", 5, "WEEK", "", "GPS week"),
    WriterField(
        "gpssec", "gpssec", ("gps", "gps_ws", "seconds"), float, "%11.3f", 11, "GPSSEC", "second", "GPS seconds"
    ),
    WriterField("apriori", "apriori", (), float, "%23.15e", 23, "APRIORI", "meter", "Apriori value parameter"),
    WriterField("estimate", "estimate", (), float, "%23.15e", 23, "ESTIMATE", "meter", "Estimated parameter value"),
    WriterField(
        "sigma", "sigma", (), float, "%11.3e", 11, "SIGMA", "meter", "Standard deviation of estimated parameter value"
    ),
    WriterField("empty", "empty", (), str, "%2s", 2, "", "", ""),  # get spaces between 'sigma' and 'param_name' column
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
        
    # Add states to WriterField depending on used pipeline
    fields_def = list(FIELDS)
    
    if dset.vars["pipeline"] == "gnss":
        
        fields_def.append(WriterField(
            "param_name",
            "param_name",
            (),
            object,
            "%-20s",
            20,
            "PARAM_NAME",
            "",
            f"Parameter name: \n"
            f"""
{'': >38}gnss_rcv_clock   - GNSS receiver clock
{'': >38}gnss_site_pos-x  - X-coordinate of site position
{'': >38}gnss_site_pos-y  - Y-coordinate of site position
{'': >38}gnss_site_pos-z  - Z-coordinate of site position
""",
        ))
        
    elif dset.vars["pipeline"] == "gnss_vel":
        
        fields_def.append(WriterField(
            "param_name",
            "param_name",
            (),
            object,
            "%-20s",
            20,
            "PARAM_NAME",
            "",
            f"Parameter name: \n"
            f"""
{'': >38}gnss_rcv_clock_drift   - GNSS receiver clock drift
{'': >38}gnss_site_vel-x        - X-coordinate of site velocity
{'': >38}gnss_site_vel-y        - Y-coordinate of site velocity
{'': >38}gnss_site_vel-z        - Z-coordinate of site velocity
""",
        ))
                
    else:
        log.fatal("Estimate solution writer is implemented only for 'gnss' and 'gnss_vel' pipeline.")
        
    # Epochwise estimation or over whole time period
    if config.tech.estimate_epochwise.bool:

        output_array = np.array([])
        for epoch in sorted(set(dset.time.gps.mjd)):
            idx = dset.time.gps.mjd == epoch

            # Append current epoch solution to final output solution for each estimated parameter
            epoch_array = _get_epoch(dset, idx, fields_def)
            output_array = np.concatenate((output_array, epoch_array), axis=0) if output_array.size else epoch_array
    else:
        # Get solution for first observation
        idx = np.squeeze(np.array(np.nonzero(dset.time.gps.mjd)) == 0)  # first observation -> TODO: Better solution?
        output_array = _get_epoch(dset, idx, fields_def)

    # Write to disk
    header = get_header(
        FIELDS,
        pgm_version=f"where {where.__version__}",
        run_by=util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else "",
        summary="Estimate solutions results",
    )
    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in fields_def),
        header=header,
        delimiter="",
        encoding="utf8",
    )


def _get_epoch(dset: "Dataset", idx: np.ndarray, fields_def: List[Tuple["WriterField"]]) -> np.ndarray:
    """Get estimation solution for given epoch indices

    Args:
        dset:       A dataset containing the data.
        idx:        Indices for epoch.
        fields_def: Writer field definition

    Returns:
        Estimation solution for given epoch indices
    """
    epoch_array = np.array([])

    # Put together fields in an array as specified by the 'dtype' tuple list
    dtype = [(f.field, f.dtype) for f in fields_def]

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
