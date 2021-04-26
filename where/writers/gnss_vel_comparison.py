"""Compare different GNSS site velocity Where datasets and write a text file

Description:
------------

A dictionary with datasets is used as input for this writer. The keys of the dictionary are station names. 

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
from midgard.data import position
from midgard.dev import plugins
from midgard.writers._writers import get_header

# Where imports
import where
from where.lib import config
from where.lib import log
from where.lib import util

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
# # PGM: where 1.1.0/midgard 1.1.2  RUN_BY:   DATE: 20210204 143718 UTC
# # DESCRIPTION: GNSS site velocity comparison results
# # 
# # 
# # 
# # HEADER         UNIT                   DESCRIPTION
# # _____________________________________________________________________________________________________________________
# # DATE                                  Date in format yyyy-mm-dd or mmm-yyyy
# # STAT                                  Station name
# # HV             m/s                    Horizontal site velocity
# # 3D             m/s                    3D site velocity
# # 
# #       DATE STAT         HV         3D
# #                        m/s        m/s
# # _____________________________________
# # 
#   2020-10-01 krss     0.0215     0.0369
#   2020-10-02 krss     0.0216     0.0347
#
# ----+----1----+----2----+----3----+----4----
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
        "station",
        object,
        "%5s",
        5,
        "STAT",
        "",
        "Station name",
    ),
    WriterField(
        "site_vel_east",
        float,
        "%11.4f",
        11,
        "EAST",
        "m/s",
        "East component of site velocity",
    ),
    WriterField(
        "site_vel_north",
        float,
        "%11.4f",
        11,
        "NORTH",
        "m/s",
        "North component of site velocity",
    ),
    WriterField(
        "site_vel_up",
        float,
        "%11.4f",
        11,
        "UP",
        "m/s",
        "Up component of site velocity",
    ),
    WriterField(
        "site_vel_h",
        float,
        "%11.4f",
        11,
        "HV",
        "m/s",
        "Horizontal site velocity",
    ),
    WriterField(
        "site_vel_3d",
        float,
        "%11.4f",
        11,
        "3D",
        "m/s",
        "3D site velocity",
    ),
)



# # PGM: where 1.1.0/midgard 1.1.2  RUN_BY:   DATE: 20210204 143718 UTC
# # DESCRIPTION: Summary of GNSS site velocity comparison results
# # 
# # HEADER         UNIT                   DESCRIPTION
# # _____________________________________________________________________________________________________________________
# # STAT                                  Station name
# # COLUMN                                Column name for which statistical analysis is carried out (e.g. site_vel_3d)
# # MEAN           m/s                    Mean for given column
# # MIN            m/s                    Minimal value for given column
# # MAX            m/s                    Maximal value for given column
# # 
# #  STAT       COLUMN     MEAN      MIN      MAX
# #                         m/s      m/s      m/s
# # _____________________________________________
# # 
#    krss   site_vel_h   0.0210   0.0192   0.0243
#    krss  site_vel_3d   0.0359   0.0330   0.0405
#    vegs   site_vel_h   0.0167   0.0158   0.0192
#    vegs  site_vel_3d   0.0318   0.0294   0.0357
#
# ----+----1----+----2----+----3----+----4----
#
FIELDS_SUM = (
    WriterField(
        "station",
        object,
        "%7s",
        5,
        "STAT",
        "",
        "Station name",
    ),
    WriterField(
        "column",
        float,
        "%16s",
        16,
        "COLUMN",
        "",
        "Column name for which statistical analysis is carried out (e.g. site_vel_3d)",
    ),
    WriterField(
        "mean",
        float,
        "%9.4f",
        9,
        "MEAN",
        "m/s",
        "Mean for given column",
    ),
    WriterField(
        "min",
        float,
        "%9.4f",
        9,
        "MIN",
        "m/s",
        "Minimal value for given column",
    ),
    WriterField(
        "max",
        float,
        "%9.4f",
        9,
        "MAX",
        "m/s",
        "Maximal value for given column",
    ),
    WriterField(
        "std",
        float,
        "%9.4f",
        9,
        "STD",
        "m/s",
        "Standard deviation for given column",
    ),
    #WriterField(
    #    "percentile",
    #    float,
    #    "%9.4f",
    #    9,
    #    "95TH",
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
def gnss_comparison(dset: "Dataset") -> None:
    """Compare GNSS datasets

    Args:
        dset:  Dictionary with station name as keys and the belonging Dataset as value
    """
    output_defs = dict()
    dset_first = dset[list(dset.keys())[0]]
    file_vars = {**dset_first.vars, **dset_first.analysis}
    file_vars["solution"] = config.tech.gnss_vel_comparison_report.solution.str.lower()
    
     # Get dataframes for writing
    _, df_day, df_month, _, _ = _generate_dataframes(dset)

    # Prepare dataframes for writing
    df_day.index = df_day.index.strftime('%Y-%m-%d')
    df_day.index.name = "date"
    df_day = df_day.reset_index()
    df_month = df_month.reset_index()
 
    # Write files for daily and monthly solutions
    samples = config.tech.gnss_vel_comparison.samples.list
    if "daily" in samples:
        output_defs.update({
              "day": df_day,
              "day_summary": _generate_dataframe_summary(df_day, index="station"),
        })

    if "monthly" in samples:
        output_defs.update({
              "month": df_month, 
              "month_summary": _generate_dataframe_summary(df_month, index="station"), 
        })

    
    for type_, output_array in output_defs.items():
        
        file_vars_tmp = file_vars.copy()
        file_vars_tmp.update(solution=f"{file_vars_tmp['solution']}_{type_}")
        file_path = config.files.path("output_gnss_vel_comparison", file_vars=file_vars_tmp)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        log.info(f"Write '{type_}' comparison file {file_path}.")
        
        fields = FIELDS_SUM if "summary" in type_ else FIELDS
        summary = "Summary of GNSS site velocity comparison results (95th percentile)" if "summary" in type_ else "GNSS site velocity comparison results (95th percentile)"
            
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
      
             
def _generate_dataframe_summary(df: pd.core.frame.DataFrame, index: str) -> pd.core.frame.DataFrame:
    """Determine dataframe summary with mean, minimal and maximal values for each numeric column
    
    Args:
        df:      Dataframe.
        index:   Chosen column of dataframe used for indexing (e.g. "station")
        
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
    """Generate dataframe based on station datasets

    The dataframe for each station in dictionary "dfs" has following columns:

        site_vel_h:   Horizontal site velocity
        site_vel_3d:  3D site velocity

    Example for "dfs" dictionary:
     
             'nabd':                       date  site_vel_h  site_vel_3d
                     0      2020-10-01 00:00:00    0.005622     0.005932
                     1      2020-10-01 00:00:00    0.005622     0.005932

             'vegs':                      date  site_vel_h  site_vel_3d
                    0      2020-10-01 00:00:00    0.005730     0.009568
                    1      2020-10-01 00:00:00    0.005730     0.009568
                      

    Example for "df_day" dictionary:

            date        site_vel_h  site_vel_3d station
            2020-10-01    0.021493     0.036855    krss
            2020-10-02    0.021583     0.034718    krss
               ...
            2020-10-30    0.018063     0.048820    nabd
            2020-10-31    0.019100     0.051127    nabd


    Example for "df_month" dictionary:

                      site_vel_h  site_vel_3d station
            date                                     
            Oct-2020    0.020967     0.035847    krss
            Oct-2020    0.016721     0.031787    vegs
            Oct-2020    0.022508     0.042471    hofs
            Oct-2020    0.020653     0.043239    hons
            Oct-2020    0.018989     0.050360    nabd


    Example for "dfs_day_field" dictionary:

            'site_vel_h':          nabf      vegs      hons      krss
                    date                                          
                    2019-03-01  1.368875  0.935687  1.136763  0.828754
                    2019-03-02  0.924839  0.728280  0.911677  0.854832


            'site_vel_3d':         nabf      vegs      hons      krss
                    date                                          
                    2019-03-01  1.715893  1.147265  1.600330  0.976541
                    2019-03-02  1.533437  1.307373  1.476295  1.136991


    Example for "dfs_month_field" dictionary:

            'site_vel_h':
                        date      nabf      vegs      hons      krss
                        Mar-2019  1.186240  0.861718  1.095827  1.021354
                        Apr-2019  0.891947  0.850343  0.977908  0.971099

            'site_vel_3d':
                        date      nabf      vegs      hons      krss
                        Mar-2019  1.854684  1.291406  1.450466  1.225467
                        Apr-2019  1.964404  1.706507  1.687994  1.500742


    Args:
        dset: Dictionary with station name as keys and the belonging Dataset as value

    Returns:
        Tuple with following entries:

        | Element              | Description                                                                          |
        |----------------------|--------------------------------------------------------------------------------------|
        | dfs                  | Dictionary with station name as keys and the belonging dataframe as value with       |
        |                      | following dataframe columns: site_vel_h, site_vel_3d                                 |
        | df_day               | Dataframe with daily entries with columns like date, stationm site_vel_3d, ...       |
        | df_month             | Dataframe with monthly entries with columns like date, stationm site_vel_3d, ...     |
        | dfs_day_field        | Dictionary with fields as keys (e.g. site_vel_3d) and the belonging dataframe as     |
        |                      | value with DAILY samples of 95th percentile and stations as columns.                 |
        | dfs_month_field      | Dictionary with fields as keys (e.g. site_vel_3d) and the belonging dataframe as     |
        |                      | value with MONTHLY samples of 95th percentile and stations as columns.               |
    """
    dsets = dset
    dfs = {}
    df_day = pd.DataFrame()
    dfs_day_field = {
        "site_vel_east": pd.DataFrame(),
        "site_vel_north": pd.DataFrame(),
        "site_vel_up": pd.DataFrame(),
        "site_vel_h": pd.DataFrame(), 
        "site_vel_3d": pd.DataFrame(), 
    }
    df_month = pd.DataFrame()
    dfs_month_field = {
        "site_vel_east": pd.DataFrame(),
        "site_vel_north": pd.DataFrame(),
        "site_vel_up": pd.DataFrame(),
        "site_vel_h": pd.DataFrame(), 
        "site_vel_3d": pd.DataFrame(), 
    }

    for station, dset in dsets.items():

        if dset.num_obs == 0:
            log.warn(f"Dataset '{station}' is empty.")
            continue

        # Determine dataframe with site_vel_h and vel_3d columns
        df = dset.as_dataframe(fields=["time.gps", "site_vel_east", "site_vel_north", "site_vel_up", "site_vel_h", "site_vel_3d"])
        df = df.rename(columns={"time_gps": "date"})
        if df.empty:
            continue
        else:
            # Save data in dictionaries
            dfs.update({station: df})

            df_day_tmp = df.set_index("date").resample("D").apply(lambda x: np.nanpercentile(x, q=95))
            for field in dfs_day_field.keys():
                if dfs_day_field[field].empty:
                    dfs_day_field[field][station] = df_day_tmp[field]
                else:
                    dfs_day_field[field] = pd.concat([dfs_day_field[field], df_day_tmp[field]], axis=1)
                dfs_day_field[field] = dfs_day_field[field].rename(columns={field: station})

            df_day_tmp["station"] = np.repeat(station, df_day_tmp.shape[0])
            df_day = pd.concat([df_day, df_day_tmp], axis=0)

            df_month_tmp = df.set_index("date").resample("M").apply(lambda x: np.nanpercentile(x, q=95))
            df_month_tmp.index = df_month_tmp.index.strftime("%b-%Y")
            for field in dfs_month_field.keys():
                dfs_month_field[field][station] = df_month_tmp[field]

            df_month_tmp["station"] = np.repeat(station, df_month_tmp.shape[0])
            df_month = pd.concat([df_month, df_month_tmp], axis=0)

    df_month.index.name = "date"

    return dfs, df_day, df_month, dfs_day_field, dfs_month_field
