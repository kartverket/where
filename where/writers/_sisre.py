"""Common functions used for writing of SISRE pipeline analysis results

Description:
------------
TODO

"""
# Standard library imports
from typing import Dict, List, Tuple

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.collections import enums

# Where imports
from where.lib import log

 
def generate_dataframes(
        dset: Dict[str, "Dataset"],
        statistics: List[str],
        samples: List[str],
        field_satellite: str, 
) -> Dict[str, pd.core.frame.DataFrame]:
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
        dset:             Dictionary with solution name as keys and the belonging Dataset as value
        statistics:       List with statistical functions, which should be applied on dataset fields. Following 
                          functions can be applied:
                              mean        - Mean
                              percentile  - 95th percentile
                              rms         - Root-mean-square
                              std         - Standard deviation
        samples:          Choose if daily and/or monthly statistical analysis dataframes should be generated.
        field_satellite:  Choose field name for which satellite based analysis should be done.
        
    Returns:
        Tuple with following entries:

        | Element              | Description                                                                          |
        |----------------------|--------------------------------------------------------------------------------------|
        | dfs                  | Dictionary with solution name as keys and the belonging dataframe as value with      |
        |                      | following dataframe columns: clk_diff, dalong_track, orb_diff_3d, sisre              |
        | dfs_day              | Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a  |
        |                      | dataframe as values with daily entries with columns like date, solution, sisre, ...  |
        | dfs_day_sat          | Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a  |
        |                      | dataframe as values with daily entries with satellite columns, which include e.g.    |
        |                      | sisre, ... solutions                                                                 |
        | dfs_month            | Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a  |
        |                      | dataframe as values with monthly entries with columns like date, solution, sisre, ...|
        | dfs_month_sat        | Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a  |
        |                      | dataframe as values with monthly entries  with satellite columns, which include e.g. |
        |                      | sisre, ... solutions                                                                 |
    """
    dsets = dset
    dfs = dict()
    dfs_day = dict()
    dfs_day_sat = dict()
    dfs_month = dict()
    dfs_month_sat = dict()
    
    fields_def = [
            "time.gps", 
            "system",
            "satellite",
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

    statistics_def = ["mean", "percentile", "rms", "std"]

    for statistic in statistics:
        if statistic not in statistics_def:
            log.fatal(f"Option '{statistic}' is not defined in 'statistic' option.")
        dfs_day.update({ statistic:  pd.DataFrame()})
        dfs_day_sat.update({ statistic:  pd.DataFrame()})
        dfs_month.update({ statistic:  pd.DataFrame()})
        dfs_month_sat.update({ statistic:  pd.DataFrame()})

    for solution, dset in dsets.items():

        if dset.num_obs == 0:
            log.warn(f"Dataset '{solution}' is empty.")
            continue
        


        # Determine dataframe 
        df = dset.as_dataframe(fields=fields_def)
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
                
                if "daily" in samples:
                    df_day_tmp, df_day_sat_tmp = _apply(df, "D", type_, field_satellite)                    
                    df_day_tmp.index = df_day_tmp.index.strftime("%Y-%m-%d")
                    _add_columns(solution, dset, df_day_tmp)    
                    dfs_day[type_] = pd.concat([dfs_day[type_], df_day_tmp], axis=0)
                    dfs_day[type_].index.name = "date"
                    
                    df_day_sat_tmp.index = df_day_sat_tmp.index.strftime("%Y-%m-%d")
                    _add_columns(solution, dset, df_day_sat_tmp) 
                    dfs_day_sat[type_] = pd.concat([dfs_day_sat[type_], df_day_sat_tmp], axis=0)
                    dfs_day_sat[type_].index.name = "date"

                if "monthly" in samples:
                    df_month_tmp, df_month_sat_tmp = _apply(df, "M", type_, field_satellite)
                    df_month_tmp.index = df_month_tmp.index.strftime("%b-%Y")
                    _add_columns(solution, dset, df_month_tmp) 
                    dfs_month[type_] = pd.concat([dfs_month[type_], df_month_tmp], axis=0)
                    dfs_month[type_].index.name = "date"
                    
                    df_month_sat_tmp.index = df_month_sat_tmp.index.strftime("%b-%Y")
                    _add_columns(solution, dset, df_month_sat_tmp) 
                    dfs_month_sat[type_] = pd.concat([dfs_month_sat[type_], df_month_sat_tmp], axis=0)
                    dfs_month_sat[type_].index.name = "date"
                      
    return dfs, dfs_day, dfs_month, dfs_day_sat, dfs_month_sat


def generate_dataframe_summary(
                    df: pd.core.frame.DataFrame, 
                    index: str, 
) -> pd.core.frame.DataFrame:
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


#
# AUXILIARY FUNCTIONS
#
def _add_columns(solution: str, dset: "Dataset", df: pd.core.frame.DataFrame) -> None:
    """Add additional columns to dataframe
    
    Args:
        solution: Solution name
        dset:     A dataset containing the data.
        df:       Dataframe. 
    """
    num_obs = df.shape[0]
    
    # Get fields to add
    frequency = None if len(dset.meta["frequencies"].keys()) > 1 else list(dset.meta["frequencies"].values())[0].replace("_", "/")
    navigation_type = None if len(dset.meta["navigation_message_type"].keys()) > 1 else list(dset.meta["navigation_message_type"].values())[0]
    service = dset.meta["service"] if "service" in dset.meta.keys() else None
    system = None if len(dset.unique("system")) > 1 else enums.gnss_id_to_name[dset.unique("system")[0]].value
    
    # Add fields as column to dataframe
    df["solution"] = np.repeat(solution, num_obs)
    if frequency:
        df["frequency"] = np.repeat(frequency, num_obs)
    if navigation_type:
        df["navigation_type"] = np.repeat(navigation_type, num_obs)  
    if service:
        df["service"] = np.repeat(service, num_obs)
    if system:
        df["system"] = np.repeat(system, num_obs)
        

def _apply(
        df: pd.core.frame.DataFrame, 
        sample: str, 
        func: str,
        field_satellite: str,
) -> Tuple[pd.core.frame.DataFrame]:
    """Resample dataframe and apply given function 

    Args:
        df:      Dataframe.
        sample:  Sample definition ("D": day, "M": month)
        func:    Function to be applied ("mean", "percentile_68", "percentile_90", "percentile", "rms", "std")
        field_satellite:  Choose field name for which satellite based analysis should be done.

    Returns:
        Tuple with resampled dataframes (date based and date+satellite based) by applying given function
    """    
    df_sat = df.pivot(index="time_gps", columns="satellite", values=field_satellite)
    df_sat_sampled = df_sat.resample(sample)
    df_reduced = df.drop([
            "system", 
            "satellite"], 
            axis=1,
    )  # Columns not needed for this analysis.
    df_sampled = df_reduced.set_index("time_gps").resample(sample)

    if func == "mean":
        df_sampled = df_sampled.mean()
        df_sat_sampled  = df_sat_sampled.mean()
    elif func == "percentile":
        df_sampled = df_sampled.apply(lambda x: np.nanpercentile(x, q=95))
        df_sat_sampled  = df_sat_sampled.apply(lambda x: np.nanpercentile(x, q=95)) 
        
        # Determine defined HAS performance metrics based on Galileo HAS SDD v1.0
        df["clk_diff_with_dt_mean_day"] = df.clk_diff_with_dt_mean - _get_daily_mean(df, "clk_diff_with_dt_mean")
        df_sampled["clk_diff_with_dt_mean_rms"] = _determine_percentile_of_epochwise_rms(
                        df, 
                        "clk_diff_with_dt_mean_day", 
                        sample,
        )
        df_sampled["orb_diff_3d_rms"] = _determine_percentile_of_epochwise_rms(df, "orb_diff_3d", sample)
        df_sampled["sisre_epoch_rms"] = _determine_percentile_of_epochwise_rms(df, "sisre", sample)
          
    elif func == "rms":
        df_sampled = df_sampled.apply(lambda x: np.sqrt(np.nanmean(np.square(x))))
        df_sat_sampled  = df_sat_sampled.apply(lambda x: np.sqrt(np.nanmean(np.square(x))))
    elif func == "std":
        df_sampled = df_sampled.std()
        df_sat_sampled = df_sat_sampled.std()
    else:
        log.fatal(f"Function '{func}' is not defined.")
        
    return df_sampled, df_sat_sampled  
   

def _determine_percentile_of_epochwise_rms(
        df: pd.core.frame.DataFrame, 
        field: str,
        sample: str,
) -> pd.core.frame.DataFrame:
    """Determine 95th percentile based on epochwise RMS solutions for given field and sample

    Args:
        df:      Dataframe.
        field:   Field name.
        sample:  Sample definition ("D": day, "M": month)

    Returns:
        95th percentile based on epochwise SISRE RMS solutions for given field and sample
    """

    # Generate monthly samples of 95th percentile FIELD based on epochwise FIELD RMS solutions
    #
    # NOTE: Following solutions assumes that each FIELD "column"-solution in dataframe 'df' is only given for one GNSS
    epochs = sorted(set(df["time_gps"]))
    df_tmp = pd.DataFrame(index=epochs, columns=[f"{field}_rms"])

    # Loop over observation epochs
    for epoch in epochs:
        idx = df["time_gps"] == epoch
        df_tmp.loc[epoch] = np.sqrt(np.nanmean(np.square(df[field][idx])))

    return df_tmp.resample(sample).apply(lambda x: np.nanpercentile(list(x), q=95))
   
             
def _get_daily_mean(
        df: pd.core.frame.DataFrame, 
        field: str,
) -> np.ndarray:
    """Get daily mean for given field

    Args:
        df:      Dataframe.
        field:   Field name.

    Returns:
        Array with daily mean of given field
    """
    dates = np.array([v.date() for v in df["time_gps"]])
    array_mean = np.zeros(df.shape[0])

    # Loop over observation epochs
    for date in sorted(set(dates)):
        idx = dates == date
        array_mean[idx] = df[field][idx].mean()

    return array_mean

