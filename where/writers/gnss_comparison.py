"""Compare different GNSS Where datasets and write a text file

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
# #       DATE STAT       EAST      NORTH         UP        HPE        VPE     3D POS   PDOP   HDOP   VDOP
# #                      meter      meter      meter      meter      meter      meter                     
# # ______________________________________________________________________________________________________
# # 
#   2020-07-01 krss     0.3300     0.4861     0.5581     0.6001     0.7522     0.8708   2.52   1.37   2.14
#   2020-07-01 vegs     0.2912     0.6225     0.6658     0.6739     0.9984     1.0891   2.57   1.25   2.31
#   2020-07-01 hons     0.4452     0.9165     0.8257     0.9482     1.0449     1.2662   2.70   1.25   2.46
# ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----
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
        "east",
        float,
        "%11.4f",
        11,
        "EAST",
        "meter",
        "xxth percentile of East component of position given in topocentric coordinate system vs. reference position",
    ),
    WriterField(
        "north",
        float,
        "%11.4f",
        11,
        "NORTH",
        "meter",
        "xxth percentile of North component of position given in topocentric coordinate system vs. reference position",
    ),
    WriterField(
        "up",
        float,
        "%11.4f",
        11,
        "UP",
        "meter",
        "xxth percentile of Up component of position given in topocentric coordinate system vs. reference position",
    ),
    WriterField(
        "hpe",
        float,
        "%11.4f",
        11,
        "HPE",
        "meter",
        "xxth percentile of Horizontal Position Error of station position vs. reference position",
    ),
    WriterField(
        "vpe",
        float,
        "%11.4f",
        11,
        "VPE",
        "meter",
        "XXth percentile of Vertical Position Error of station position vs. reference position",
    ),
    WriterField(
        "pos_3d",
        float,

        "%11.4f",
        11,
        "3D POS",
        "meter",
        "xxth percentile of 3D position error of station position vs. reference position",
    ),
    WriterField("pdop", float, "%7.2f", 7, "PDOP", "", "Position dilution of precision"),
    WriterField("hdop", float, "%7.2f", 7, "HDOP", "", "Horizontal dilution of precision"),
    WriterField("vdop", float, "%7.2f", 7, "VDOP", "", "Vertical dilution of precision"),
)




# # PGM: where 1.0.4/midgard 1.1.1  RUN_BY: NMA  DATE: 20201022 203257 UTC
# # DESCRIPTION: Summary of GNSS comparison results
# # 
# #  STAT  COLUMN     MEAN      MIN      MAX
# #                  meter    meter    meter
# # ________________________________________
# # 
#    krss    pdop   2.5236   2.5236   2.5236
#    krss    hdop   1.3720   1.3720   1.3720
#    krss    vdop   2.1351   2.1351   2.1351
#    krss    east   0.3300   0.3300   0.3300
#    krss   north   0.4861   0.4861   0.4861
#    krss      up   0.5581   0.5581   0.5581
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
        "%8s",
        8,
        "COLUMN",
        "",
        "Column name for which statistical analysis is carried out (e.g. HPE, VPE)",
    ),
    WriterField(
        "mean",
        float,
        "%9.4f",
        9,
        "MEAN",
        "meter",
        "Mean for given column",
    ),
    WriterField(
        "min",
        float,
        "%9.4f",
        9,
        "MIN",
        "meter",
        "Minimal value for given column",
    ),
    WriterField(
        "max",
        float,
        "%9.4f",
        9,
        "MAX",
        "meter",
        "Maximal value for given column",
    ),
    WriterField(
        "std",
        float,
        "%9.4f",
        9,
        "STD",
        "meter",
        "Standard deviation for given column",
    ),
    #WriterField(
    #    "percentile",
    #    float,
    #    "%9.4f",
    #    9,
    #    "xxth",
    #    "meter",
    #    "xxth percentile for given column",
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
                  "day_summary": _generate_dataframe_summary(dfs_day[type_], index="station"),
            })

        if "monthly" in samples:
            output_defs.update({
                  "month": dfs_month[type_], 
                  "month_summary": _generate_dataframe_summary(dfs_month[type_], index="station"), 
            })

        for sample, output_array in output_defs.items():
            
            file_vars_tmp = file_vars.copy()
            file_vars_tmp.update(solution=f"{file_vars_tmp['solution']}_{type_}_{sample}")
            file_path = config.files.path("output_gnss_comparison", file_vars=file_vars_tmp)
            file_path.parent.mkdir(parents=True, exist_ok=True)
            
            log.info(f"Write '{sample}' comparison file {file_path}.")
            
            fields = FIELDS_SUM if "summary" in sample else FIELDS
            summary = "Summary of GNSS comparison results (xxth percentile)" if "summary" in sample else "GNSS comparison results (xxth percentile)"
                
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
        func:    Function to be applied ("mean", "percentile_68", "percentile_90", "percentile_95", "rms", "std")

    Returns:
        Resampled dataframe by applying given function
    """
    df_sampled = df.set_index("time_gps").resample(sample)

    if func == "mean":
        df_sampled = df_sampled.mean()
    elif func == "percentile_68":
        df_sampled = df_sampled.apply(lambda x: np.nanpercentile(x, q=68))
    elif func == "percentile_90":
        df_sampled = df_sampled.apply(lambda x: np.nanpercentile(x, q=90))
    elif func == "percentile_95":
        df_sampled = df_sampled.apply(lambda x: np.nanpercentile(x, q=95))
    elif func == "rms":
        df_sampled = df_sampled.apply(lambda x: np.sqrt(np.nanmean(np.square(x))))
    elif func == "std":
        df_sampled = df_sampled.std()
    else:
        log.fatal(f"Function '{func}' is not defined.")

    return df_sampled  
      
             
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

        east:   East-coordinate in topocentric system
        north:  North-coordinate in topocentric system
        up:     Up-coordinate in topocentric system
        hpe:    horizontal position error
        vpe:    vertical position error
        pos_3d: 3D position error
        pdop:   position dilution of precision
        hdop:   horizontal dilution of precision
        vdop:   vertical dilution of precision

    Example for "dfs" dictionary:
     
            'hons':                   date       hpe       vpe      east     north        up
                    0      2019-03-01 00:00:00  0.301738  0.057244  0.113758  0.279472  0.057244
                    1      2019-03-01 00:00:00  0.301738  0.057244  0.113758  0.279472  0.057244

            'krss':                   date       hpe       vpe      east     north        up
                    0      2019-03-01 00:00:00  0.710014  0.186791 -0.235267  0.669903  0.186791
                    1      2019-03-01 00:00:00  0.710014  0.186791 -0.235267  0.669903  0.186791


    Example for "dfs_day" dictionary for "mean" key:
        'mean':{
                            pdop      hdop      vdop  ...       vpe    pos_3d  station
            date                                      ...                             
            2020-07-01  2.523569  1.371987  2.135124  ...  0.752227  0.870759     krss
            2020-07-01  2.571588  1.247443  2.308469  ...  0.998428  1.089063     vegs
            2020-07-01  2.622492  1.289113  2.330969  ...  1.084772  1.220454     hofs
            2020-07-01  2.699645  1.246052  2.456847  ...  1.044877  1.266227     hons
            2020-07-01  2.695779  1.156999  2.448314  ...  1.461449  1.619489     nabd
        }


    Example for "dfs_month" dictionary for "mean" key:
        'mean':{
                          pdop      hdop      vdop  ...       vpe    pos_3d  station
            date
            Jul-2020  2.523569  1.371987  2.135124  ...  0.752227  0.870759     krss
            Jul-2020  2.571588  1.247443  2.308469  ...  0.998428  1.089063     vegs
            Jul-2020  2.622492  1.289113  2.330969  ...  1.084772  1.220454     hofs
            Jul-2020  2.699645  1.246052  2.456847  ...  1.044877  1.266227     hons
            Jul-2020  2.695779  1.156999  2.448314  ...  1.461449  1.619489     nabd
        }

    Args:
        dset: Dictionary with station name as keys and the belonging Dataset as value

    Returns:
        Tuple with following entries:

        | Element              | Description                                                                          |
        |----------------------|--------------------------------------------------------------------------------------|
        | dfs                  | Dictionary with station name as keys and the belonging dataframe as value with       |
        |                      | following dataframe columns: east, north, up, hpe, vpe, pos_3d                       |
        | dfs_day              | Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a  |
        |                      | dataframe as values with daily entries with columns like date, station, hpe, vpe, ...|
        | dfs_month            | Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a  |
        |                      | dataframe as values with monthly entries with columns like date, station, hpe,       |
        |                      | vpe, ...                                                                             |
    """
    dsets = dset
    dfs = dict()
    dfs_day = dict()
    dfs_month = dict()

    statistics_def = ["mean", "percentile_68", "percentile_90", "percentile_95", "rms", "std"]
    statistics_cfg = config.tech[_SECTION].statistic.list

    for statistic in statistics_cfg:
        if statistic not in statistics_def:
            log.fatal(f"Option '{statistic}' is not defined in 'statistic' option.")
        dfs_day.update({ statistic:  pd.DataFrame()})
        dfs_month.update({ statistic:  pd.DataFrame()})

    for station, dset in dsets.items():

        if dset.num_obs == 0:
            log.warn(f"Dataset '{station}' is empty.")
            continue

        # Determine topocentric coordinates (east, north, up)
        ref_pos = position.Position(
                            np.repeat(
                                  np.array([dset.meta["pos_x"], dset.meta["pos_y"], dset.meta["pos_z"]])[None, :], 
                                  dset.num_obs, 
                                  axis=0,
                            ), 
                            system="trs",
        )

        if not "enu" in dset.fields:
            dset.add_position_delta(
                name="enu",
                val=(dset.site_pos.trs - ref_pos).val,
                system="trs",
                ref_pos=ref_pos,
                write_level="operational",
            )

        # TODO: Maybe it is not necessary to introduce enu, hpe and vpe to dataset
        #      Maybe better to introduce fields in estimate stage already.
        if not "hpe" in dset.fields:
            hpe = np.sqrt(dset.enu.enu.east ** 2 + dset.enu.enu.north ** 2)
            dset.add_float("hpe", val=hpe, write_level="operational")

        if not "vpe" in dset.fields:
            vpe = np.absolute(dset.enu.enu.up)
            dset.add_float("vpe", val=vpe, write_level="operational")

        if not "pos_3d" in dset.fields:
            pos_3d = np.sqrt(dset.enu.enu.east ** 2 + dset.enu.enu.north ** 2 + dset.enu.enu.up ** 2)
            dset.add_float("pos_3d", val=pos_3d, write_level="operational")

        # Determine dataframe 
        df = dset.as_dataframe(fields=["enu.enu", "time.gps", "hpe", "vpe", "pos_3d", "pdop", "vdop", "hdop"])
        df = df.rename(columns={"enu_enu_0": "east", "enu_enu_1": "north", "enu_enu_2": "up"})

        if df.empty:
            continue
        else:
            # Save data in dictionaries
            dfs.update({station: df})

            for type_ in dfs_day.keys():

                df_day_tmp = _apply(df, "D", type_)
                df_day_tmp.index = df_day_tmp.index.strftime("%Y-%m-%d")
                df_day_tmp["station"] = np.repeat(station, df_day_tmp.shape[0])
                dfs_day[type_] = pd.concat([dfs_day[type_], df_day_tmp], axis=0)
                dfs_day[type_].index.name = "date"

                df_month_tmp = _apply(df, "M", type_)
                df_month_tmp.index = df_month_tmp.index.strftime("%b-%Y")
                df_month_tmp["station"] = np.repeat(station, df_month_tmp.shape[0])
                dfs_month[type_] = pd.concat([dfs_month[type_], df_month_tmp], axis=0)
                dfs_month[type_].index.name = "date"

    return dfs, dfs_day, dfs_month
