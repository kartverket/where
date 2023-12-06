"""Compare different SISRE Where datasets and write a text file

Description:
------------

A dictionary with datasets is used as input for this writer. The keys of the dictionary are solution names. 

Example:
--------


"""
# Standard library imports
from collections import namedtuple
from typing import Dict

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.dev import plugins
from midgard.writers._writers import get_header

# Where imports
import where
from where.lib import config
from where.lib import log
from where.lib import util

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
        "%30s",
        30,
        "SOLUTION",
        "",
        "Solution name",
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
        "Satellite clock correction difference related to center of mass of satellite and corrected for "
        "satellite biases",
    ),
    WriterField(
        "clk_diff_with_dt_mean",
        float,
        "%16.4f",
        16,
        "ΔCLOCK_MEAN",
        "meter",
        "Satellite clock correction difference related to center of mass of satellite and corrected for "
        "satellite biases and averaged clock offset in each epoch",
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
        "%30s",
        28,
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
    file_vars = {**dset_first.vars, **dset_first.analysis}
    file_vars["solution"] = config.tech[_SECTION].solution.str.lower()
    samples = config.tech[_SECTION].samples.list
    
     # Get dataframes for writing
    _, dfs_day, dfs_month = _generate_dataframes(dset)

    # Write files for daily and monthly solutions
    for type_ in dfs_day.keys():

        # Prepare dataset for writing
        dfs_day[type_] = dfs_day[type_].reset_index()
        dfs_month[type_] = dfs_month[type_].reset_index()
        
        if "daily" in samples:
            output_defs.update({
                  "day": dfs_day[type_],
                  "day_summary": _generate_dataframe_summary(dfs_day[type_], index="solution"),
            })
        if "monthly" in samples:
            output_defs.update({
                  "month": dfs_month[type_], 
                  "month_summary": _generate_dataframe_summary(dfs_month[type_], index="solution"), 
            })

        for sample, output_array in output_defs.items():
            
            file_vars_tmp = file_vars.copy()
            file_vars_tmp.update(solution=f"{file_vars_tmp['solution']}_{type_}_{sample}")
            file_path = config.files.path("output_sisre_comparison", file_vars=file_vars_tmp)
            file_path.parent.mkdir(parents=True, exist_ok=True)
            
            log.info(f"Write '{sample}' comparison file {file_path}.")
            
            fields = FIELDS_SUM if "summary" in sample else FIELDS
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


#
# AUXILIARY FUNCTIONS
#
def _apply(df: pd.core.frame.DataFrame, sample: str, func: str) -> pd.core.frame.DataFrame:
    """Resample dataframe and apply given function 

    Args:
        df:      Dataframe.
        sample:  Sample definition ("D": day, "M": month)
        func:    Function to be applied ("mean", "percentile_68", "percentile_90", "percentile", "rms", "std")

    Returns:
        Resampled dataframe by applying given function
    """
    df_sampled = df.set_index("time_gps").resample(sample)

    if func == "mean":
        df_sampled = df_sampled.mean()
    elif func == "percentile":
        df_sampled = df_sampled.apply(lambda x: np.nanpercentile(x, q=95))
        df_sampled["sisre_epoch_rms"] = _determine_sisre_percentile_of_epochwise_rms(df, sample)
        
    elif func == "rms":
        df_sampled = df_sampled.apply(lambda x: np.sqrt(np.nanmean(np.square(x))))
    elif func == "std":
        df_sampled = df_sampled.std()
    else:
        log.fatal(f"Function '{func}' is not defined.")

    return df_sampled  
   

def _determine_sisre_percentile_of_epochwise_rms(df: pd.core.frame.DataFrame, sample: str) -> pd.core.frame.DataFrame:
    """Determine 95th percentile of SISRE based on epochwise SISRE RMS solutions for given sample

    Args:
        df:      Dataframe.
        sample:  Sample definition ("D": day, "M": month)

    Returns:
        95th percentile of SISRE based on epochwise SISRE RMS solutions for given sample
    """

    # Generate monthly samples of 95th percentile FIELD based on epochwise FIELD RMS solutions(after SDD v1.1 version)
    #
    # NOTE: Following solutions assumes that each FIELD "column"-solution in dataframe 'df' is only given for one GNSS
    epochs = sorted(set(df["time_gps"]))
    df_tmp = pd.DataFrame(index=epochs, columns=["sisre_epoch_rms"])

    # Loop over observation epochs
    for epoch in epochs:
        idx = df["time_gps"] == epoch
        df_tmp.loc[epoch] = np.sqrt(np.nanmean(np.square(df.sisre[idx])))

    return df_tmp.resample(sample).apply(lambda x: np.nanpercentile(list(x), q=95))
   
             
def _generate_dataframe_summary(df: pd.core.frame.DataFrame, index: str) -> pd.core.frame.DataFrame:
    """Determine dataframe summary with mean, minimal and maximal values for each numeric column
    
    Args:
        df:      Dataframe.
        index:   Chosen column of dataframe used for indexing (e.g. "solution")
        
    Returns:
        Dataframe with statistics (mean, min, max, std, percentile, rms) based on given numeric dataframe columns
    """
    functions = ["mean", "min", "max", "std", "percentile", "rms"]
    ncols = df.select_dtypes("number").columns  # only use of numeric columns
    items = df[index].unique()
    
    # Initialize summary dataframe with chosen 'index' column and column names
    summary = pd.DataFrame(columns=[index, "column"] + functions)
    entries = []
    columns= []
    for item in items: 
        entries.extend([item] * len(ncols)) 
        columns.extend(ncols)
    summary[index] = entries
    summary["column"] = columns
    
    # Determine statistics 
    for item in items:
        idx_df = df[index] == item
        idx_sum = summary[index] == item
        
        for func in functions:
            if func == "rms":
                summary[func][idx_sum] = df[ncols][idx_df].apply(lambda x: np.sqrt(np.nanmean(np.square(x))))
            elif func == "percentile":
                summary[func][idx_sum] = df[ncols][idx_df].apply(lambda x: np.nanpercentile(x, q=95))
            else:
                summary[func][idx_sum] = getattr(df[ncols][idx_df], func)().values  # mean, max, min and std dataframe
                                                                                    # functions skip NaN values
            
    return summary
              
    
def _generate_dataframes(dset: Dict[str, "Dataset"]) -> Dict[str, pd.core.frame.DataFrame]:
    """Generate dataframe based on solution datasets

    The dataframe for each solution in dictionary "dfs" has following columns:

        age_of_ephemeris:        Age of ephemeris, which is the difference between the observation time and the time
                                 of ephemeris (ToE)
        clk_diff:                Satellite clock correction difference related to center of mass of satellite and
                                 corrected for satellite biases
        clk_diff_with_dt_mean:   Satellite clock correction difference related to center of mass of satellite and 
                                 corrected for satellite biases and averaged clock offset in each epoch       
        dalong_track:            Satellite coordinate difference between broadcast and precise ephemeris in
                                 along-track direction
        dcross_track:            Satellite coordinate difference between broadcast and precise ephemeris in
                                 cross-track direction
        dradial:                 Satellite coordinate difference between broadcast and precise ephemeris in
                                 radial direction
        orb_diff_3d:             3D orbit difference
        sisre_orb:               Orbit-only SISRE
        sisre:                   Global averaged SISRE  
     
    Example for "dfs" dictionary:
     
            '_gal_has_concatenated':               time_gps     clk_diff  dalong_track  dcross_track dradial   sisre
                                       0   2019-03-01 00:00:00  0.301738  0.057244      0.113758     0.279472  0.057244
                                       1   2019-03-01 00:00:00  0.301738  0.057244      0.113758     0.279472  0.057244

            '_gal_os_concatenated':                time_gps     clk_diff  dalong_track  dcross_track dradial   sisre
                                       0   2019-03-01 00:00:00  0.710014  0.186791 -0.235267  0.669903  0.186791
                                       1   2019-03-01 00:00:00  0.710014  0.186791 -0.235267  0.669903  0.186791


    Example for "dfs_day" dictionary for "mean" key:
        'mean':{
                        clk_diff  dalong_track dcross_track ...  dradial   sisre  solution
            date                                      ...                             
            2020-07-01  2.523569  1.371987  2.135124  ...  0.752227  0.870759     _gal_has_concatenated
            2020-07-01  2.571588  1.247443  2.308469  ...  0.998428  1.089063     _gal_os_concatenated
            2020-07-01  2.622492  1.289113  2.330969  ...  1.084772  1.220454     _gps_has_concatenated
            2020-07-01  2.699645  1.246052  2.456847  ...  1.044877  1.266227     _gps_os_concatenated

        }


    Example for "dfs_month" dictionary for "mean" key:
                        clk_diff  dalong_track dcross_track ...  dradial   sisre  solution
            date                                      ...                             
            Jul-2020    2.523569  1.371987  2.135124  ...  0.752227  0.870759     _gal_has_concatenated
            Jul-2020    2.571588  1.247443  2.308469  ...  0.998428  1.089063     _gal_os_concatenated
            Jul-2020    2.622492  1.289113  2.330969  ...  1.084772  1.220454     _gps_has_concatenated
            Jul-2020    2.699645  1.246052  2.456847  ...  1.044877  1.266227     _gps_os_concatenated
        }

    Args:
        dset: Dictionary with solution name as keys and the belonging Dataset as value

    Returns:
        Tuple with following entries:

        | Element              | Description                                                                          |
        |----------------------|--------------------------------------------------------------------------------------|
        | dfs                  | Dictionary with solution name as keys and the belonging dataframe as value with      |
        |                      | following dataframe columns: clk_diff, dalong_track, orb_diff_3d, sisre              |
        | dfs_day              | Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a  |
        |                      | dataframe as values with daily entries with columns like date, solution, sisre, ...  |
        | dfs_month            | Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a  |
        |                      | dataframe as values with monthly entries with columns like date, solution, sisre, ...|
    """
    dsets = dset
    dfs = dict()
    dfs_day = dict()
    dfs_month = dict()

    statistics_def = ["mean", "percentile", "rms", "std"]
    statistics_cfg = config.tech[_SECTION].statistic.list

    for statistic in statistics_cfg:
        if statistic not in statistics_def:
            log.fatal(f"Option '{statistic}' is not defined in 'statistic' option.")
        dfs_day.update({ statistic:  pd.DataFrame()})
        dfs_month.update({ statistic:  pd.DataFrame()})

    for solution, dset in dsets.items():

        if dset.num_obs == 0:
            log.warn(f"Dataset '{solution}' is empty.")
            continue

        # Determine dataframe 
        df = dset.as_dataframe(fields=[
                    "time.gps", 
                    "age_of_ephemeris",
                    "clk_diff", 
                    "clk_diff_with_dt_mean", 
                    "orb_diff.acr.along",
                    "orb_diff.acr.cross",
                    "orb_diff.acr.radial",
                    "orb_diff_3d", 
                    "sisre_orb", 
                    "sisre",
            ]
        )
        df = df.rename(columns={
                    "orb_diff_acr_along": "dalong_track", 
                    "orb_diff_acr_cross": "dcross_track", 
                    "orb_diff_acr_radial": "dradial",
        })

        if df.empty:
            continue
        else:
            # Save data in dictionaries
            dfs.update({solution: df})

            for type_ in dfs_day.keys():

                df_day_tmp = _apply(df, "D", type_)
                df_day_tmp.index = df_day_tmp.index.strftime("%Y-%m-%d")
                df_day_tmp["solution"] = np.repeat(solution, df_day_tmp.shape[0])
                dfs_day[type_] = pd.concat([dfs_day[type_], df_day_tmp], axis=0)
                dfs_day[type_].index.name = "date"

                df_month_tmp = _apply(df, "M", type_)
                df_month_tmp.index = df_month_tmp.index.strftime("%b-%Y")
                df_month_tmp["solution"] = np.repeat(solution, df_month_tmp.shape[0])
                dfs_month[type_] = pd.concat([dfs_month[type_], df_month_tmp], axis=0)
                dfs_month[type_].index.name = "date"
  
    return dfs, dfs_day, dfs_month
