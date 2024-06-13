"""Compare different SISRE Where datasets and write a text file

Description:
------------

A dictionary with datasets is used as input for this writer. The keys of the dictionary are solution names. 

Example:
--------


"""
# Standard library imports
from collections import namedtuple, OrderedDict
from typing import List, Tuple

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.dev import plugins
from midgard.writers._writers import get_header
from midgard.writers.csv_ import csv_

# Where imports
import where
from where.lib import config, log, util
from where.writers._sisre import generate_dataframes, generate_dataframe_summary

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])

WriterField = namedtuple(
    "WriterField", ["name", "dtype", "format", "width", "header", "unit", "description"]
)
WriterField.__new__.__defaults__ = (None,) * len(WriterField._fields)
WriterField.__doc__ = """A convenience class for defining a output field for the writer

    Args:
        name  (str):             Field name
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
# #       DATE                      SOLUTION             AGE          ΔCLOCK     ΔCLOCK_MEAN    ΔALONG_TRACK   ...
# #                                                   second           meter           meter           meter  
# # _______________________________________________________________________________________________________________
#   2023-05-01         _gal_has_concatenated            4200          0.2763          0.1441          0.1465   ...
#   2023-05-01          _gal_os_concatenated            3900          0.3424          0.2572          0.3692   ...
# ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1---
#
FIELDS = (
    WriterField(
        "date",
        object,
        "%12s",
        10,
        "DATE",
        "",
        "Date in format yyyy-mm-dd or mmm-yyyy",
    ),
    WriterField(
        "solution",
        object,
        "%40s",
        40,
        "SOLUTION",
        "",
        "Solution name",
    ),
    WriterField(
        "system",
        object,
        "%10s",
        10,
        "SYS",
        "",
        "GNSS system",
    ),
    WriterField(
        "service",
        object,
        "%6s",
        6,
        "SRVC",
        "",
        "GNSS service",
    ),
    WriterField(
        "satellite",
        object,
        "%5s",
        5,
        "SAT",
        "",
        "Satellite identifier",
    ),
    WriterField(
        "navigation_type",
        object,
        "%6s",
        6,
        "NAV",
        "",
        "Navigation message type used in analysis",
    ),
    WriterField(
        "frequency",
        object,
        "%8s",
        8,
        "FREQ",
        "",
        "Frequency used in analysis",
    ),
    WriterField(
        "age_of_ephemeris",
        float,
        "%16.0f",
        16,
        "AGE",
        "second",
        "Age of ephemeris, which is the difference between the observation time and the time of ephemeris (ToE)",
    ),
    WriterField(
        "clk_diff",
        float,
        "%16.4f",
        16,
        "ΔCLOCK",
        "meter",
        "Satellite clock correction difference corrected for satellite biases",
    ),
    WriterField(
        "clk_diff_with_dt_mean",
        float,
        "%16.4f",
        16,
        "ΔCLOCK_MEAN",
        "meter",
        "Satellite clock correction difference corrected for satellite biases and averaged clock offset in each epoch",
    ),
    WriterField(
        "clk_diff_with_dt_mean_rms",
        float,
        "%16.4f",
        16,
        "ΔCLOCK_MEAN_RMS",
        "meter",
        "Satellite clock correction difference corrected for satellite biases and averaged clock offset in each epoch. "
        "Afterwards the epochwise RMS is determined.",
    ),
    WriterField(
        "clk_diff_with_dt_mean_day_rms",
        float,
        "%16.4f",
        16,
        "ΔCLOCK_MDAY_RMS",
        "meter",
        "Satellite clock correction difference corrected for satellite biases and averaged clock offset in each epoch. "
        "In addition RMS for each satellite and daily samples is subtracted and afterwards the epochwise RMS is "
        "determined.",
    ),
    WriterField(
        "dalong_track",
        float,
        "%16.4f",
        16,
        "ΔALONG_TRACK",
        "meter",
        "Satellite coordinate difference between broadcast and precise ephemeris in along-track direction",
    ),
    WriterField(
        "dcross_track",
        float,
        "%16.4f",
        16,
        "ΔCROSS_TRACK",
        "meter",
        "Satellite coordinate difference between broadcast and precise ephemeris in cross-track direction",
    ),
    WriterField(
        "dradial",
        float,
        "%16.4f",
        16,
        "ΔRADIAL",
        "meter",
        "Satellite coordinate difference between broadcast and precise ephemeris in radial direction",
    ),
    WriterField("orb_diff_3d", float, "%16.4f", 16, "ORB_DIFF_3D", "meter", "3D orbit difference"),
    WriterField(
        "orb_diff_3d_rms", 
        float, 
        "%16.4f", 
        16, 
        "ORB_DIFF_3D_RMS", 
        "meter", 
        "3D orbit difference based on epochwise RMS solution",
    ),
    WriterField("sisre_orb", float, "%16.4f", 16, "SISRE_ORB", "meter", "Orbit-only SISRE"),
    WriterField(
        "sisre", 
        float, 
        "%16.4f", 
        16, 
        "SISRE", 
        "meter", 
        "Global averaged SISRE",
    ),
    WriterField(
        "sisre_epoch_rms", 
        float, 
        "%16.4f", 
        16, 
        "SISRE_RMS", 
        "meter", 
        "Global averaged SISRE based on epochwise RMS solutions",
    ),
)




# # PGM: where 2.1.3/midgard 1.2.2  RUN_BY: NMA  DATE: 20231206 154952 UTC
# # DESCRIPTION: Summary of SISRE comparison results (95th percentile)
# # 
# #                     SOLUTION                        COLUMN       MEAN        MIN        MAX        STD
# #                                                                 meter      meter      meter      meter
# # ______________________________________________________________________________________________________
# # 
#          _gal_has_concatenated                      clk_diff     0.2333     0.1902     0.2763     0.0609
#          _gal_has_concatenated         clk_diff_with_dt_mean     0.1480     0.1441     0.1519     0.0055
#          _gal_has_concatenated                         sisre     0.1685     0.1629     0.1742     0.0080
# ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----
#
FIELDS_SUM = (
    WriterField(
        "solution",
        object,
        "%40s",
        38,
        "SOLUTION",
        "",
        "Solution name",
    ),
    WriterField(
        "column",
        float,
        "%30s",
        30,
        "COLUMN",
        "",
        "Column name for which statistical analysis is carried out (e.g. SISRE)",
    ),
    WriterField(
        "mean",
        float,
        "%11.4f",
        11,
        "MEAN",
        "meter",
        "Mean for given column",
    ),
    WriterField(
        "min",
        float,
        "%11.4f",
        11,
        "MIN",
        "meter",
        "Minimal value for given column",
    ),
    WriterField(
        "max",
        float,
        "%11.4f",
        11,
        "MAX",
        "meter",
        "Maximal value for given column",
    ),
    WriterField(
        "std",
        float,
        "%11.4f",
        11,
        "STD",
        "meter",
        "Standard deviation for given column",
    ),
    #WriterField(
    #    "percentile",
    #    float,
    #    "%9.4f",
    #    9,
    #    "95th",
    #    "meter",
    #    "95th percentile for given column",
    #),
    #WriterField(
    #    "rms",
    #    float,
    #    "%9.4f",
    #    9,
    #    "RMS",
    #    "meter",
    #    "Root mean square for given column",
    #),
)


@plugins.register
def sisre_comparison(dset: "Dataset") -> None:
    """Compare SISRE datasets

    Args:
        dset:  Dictionary with solution name as keys and the belonging Dataset as value
    """
    output_defs = dict()
    dset_first = dset[list(dset.keys())[0]]
    samples = config.tech[_SECTION].samples.list

    #+WORKAROUND_MURKS: service meta information is not available in earlier datasets.
    for sol in dset.keys():
        if "service" not in dset[sol].meta.keys():
            service = "has" if "has" in sol else "os"
            dset[sol].meta["service"] = service.upper()
    #-WORKAROUND_MURKS:

    # Add service information (OS and/or HAS) to file variables
    file_vars = {**dset_first.vars, **dset_first.analysis}
    for sol in dset.keys(): 

        if "HAS" in dset[sol].meta["service"]:
            file_vars["service"] = "has"
            file_vars["SERVICE"] = "HAS"
            break
        else:
            file_vars["service"] = "os"
            file_vars["SERVICE"] = "OS"

    # CSV file options
    mode_csv = config.tech[_SECTION].mode_csv.str
    write_csv = config.tech[_SECTION].write_csv.bool
            
     # Get dataframes for writing
    _, dfs_day, dfs_month, dfs_day_sat, dfs_month_sat = generate_dataframes(
                                                        dset,
                                                        statistics=config.tech[_SECTION].statistic.list,
                                                        samples=samples,
                                                        field_satellite=config.tech[_SECTION].field_satellite.str,
    )
    
    # Write files for daily and monthly solutions by looping over statistical solutions (e.g. type_='percentile')
    for type_ in dfs_day.keys():

        # Prepare dataset for writing
        dfs_day[type_] = dfs_day[type_].reset_index()
        dfs_day_sat[type_] = dfs_day_sat[type_].reset_index()
        dfs_month[type_] = dfs_month[type_].reset_index()
        dfs_month_sat[type_] = dfs_month_sat[type_].reset_index()
        
        if "daily" in samples:
            if dfs_day[type_].empty:
                log.warn(f"Daily '{type_}' dataframe is empty.")
                
            else:
                output_defs.update({
                      "day": dfs_day[type_],
                      "day_sat": dfs_day_sat[type_],
                      "day_summary": generate_dataframe_summary(dfs_day[type_], index="solution"),
                })
        if "monthly" in samples:
            if dfs_month[type_].empty:
                log.warn(f"Monthly '{type_}' dataframe is empty.")
            
            else:
                output_defs.update({
                      "month": dfs_month[type_], 
                      "month_sat": dfs_month_sat[type_],
                      "month_summary": generate_dataframe_summary(dfs_month[type_], index="solution"), 
                })
                           
        for sample, output_array in output_defs.items():
            
            file_vars_tmp = file_vars.copy()
            file_vars_tmp.update(solution=f"{type_}_{sample}")
            file_path = config.files.path("output_sisre_comparison", file_vars=file_vars_tmp)
            file_path.parent.mkdir(parents=True, exist_ok=True)
            
            log.info(f"Write '{sample}' comparison file {file_path}.")
            
            fields = FIELDS_SUM if "summary" in sample else FIELDS
            fields = _get_existing_fields(output_array.columns, fields) 
            summary = "Summary of SISRE comparison results (95th percentile)" if "summary" in sample else "SISRE comparison results (95th percentile)"
                            
            # Get header
            header = get_header(
                fields,
                pgm_version=f"where {where.__version__}",
                run_by=util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else "",
                summary=summary,
            )

            # Write to disk
            np.savetxt(
                file_path,
                output_array[[f.name for f in fields]].to_numpy(),
                fmt=tuple(f.format for f in fields),
                header=header,
                delimiter="",
                encoding="utf8",
            )
            
            # Write daily and monthly CSV files
            if write_csv:
                file_vars_tmp = file_vars.copy()
                file_vars_tmp.update(solution=f"{type_}_{sample}")
                file_path = config.files.path("output_sisre_comparison_csv", file_vars=file_vars_tmp)
                file_path.parent.mkdir(parents=True, exist_ok=True)
                log.info(f"Write '{sample}' comparison file {file_path} in CSV format.")
                output_array.to_csv(file_path, float_format="%.3f", index=False, na_rep="nan", mode=mode_csv)


#
# AUXILIARY FUNCTIONS
#
def _get_existing_fields(data_fields: List[str], writer_fields: Tuple["WriterField", ...]) -> List["WriterField"]:
    """Get existing writer fields, which are given in Dataset.

    Args:
        data_fields:      List with data fields
        writers_fields:   Writer fields

    Returns:
        Existing writer fields
    """
    common_fields = []
    for writer in writer_fields:
        if writer.name in data_fields:
            common_fields.append(writer)
    return common_fields


def _get_sat_fields(df: pd.core.frame.DataFrame) -> Tuple[WriterField]:
    """Get satellite WriterField definition
    
    Args:
        df: Dataframe
        
    Returns:
        Tuple with WriterField definition for each satellite column
    """
    satellites = set(df.columns) - set(["date", "frequency", "navigation_type", "service", "solution", "system"])
    field_sat = ([
        WriterField(
            "date",
            object,
            "%12s",
            10,
            "DATE",
            "",
            "Date in format yyyy-mm-dd or mmm-yyyy",
        ),
        WriterField(
            "solution",
            object,
            "%30s",
            28,
            "SOLUTION",
            "",
            "Solution name",
        ),
        WriterField(
            "system",
            object,
            "%10s",
            10,
            "SYS",
            "",
            "GNSS system",
        ),
        WriterField(
            "service",
            object,
            "%6s",
            6,
            "SRVC",
            "",
            "GNSS service",
        ),
        WriterField(
            "navigation_type",
            object,
            "%6s",
            6,
            "NAV",
            "",
            "Navigation message type used in analysis",
        ),
        WriterField(
            "frequency",
            object,
            "%8s",
            8,
            "FREQ",
            "",
            "Frequency used in analysis",
        ),
    ])
       
    for sat in sorted(satellites):
        field_sat = field_sat + ([
            WriterField(
                sat,
                float,
                "%11.4f",
                11,
                sat,
                "meter",
                f"95th percentile SISRE for satellite {sat}",
            )
    ])
        
    return field_sat
