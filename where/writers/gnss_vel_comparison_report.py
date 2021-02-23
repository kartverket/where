"""Compare different GNSS site velocity Where datasets

Description:
------------

A dictionary with datasets is used as input for this writer. The keys of the dictionary are station names. 

Example:
--------

    from where import data
    from where import writers

    # Read a dataset
    dset = data.Dataset(rundate=rundate, tech=tech, stage=stage, dataset_name=name, dataset_id=dataset_id)

    # Write dataset
    writers.write_one('gnss_vel_comparison_report', dset=dset, do_report=False)

"""
# Standard library imports
import copy
from typing import Any, Dict
from pathlib import PosixPath

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.dev import plugins
from midgard.plot.matplotlib_extension import plot

# Where imports
from where.lib import config
from where.lib import log
from where.postprocessors.gnss_velocity_fields import gnss_velocity_fields
from where.writers._report import Report

FIGURE_FORMAT = "png"
FILE_NAME = __name__.split(".")[-1]


@plugins.register
def gnss_vel_comparison_report(dset: Dict[str, "Dataset"]) -> None:
    """Compare GNSS site velocity datasets

    Args:
        dset:  Dictionary with station name as keys and the belonging Dataset as value
    """
    dset_first = dset[list(dset.keys())[0]]
    file_vars = {**dset_first.vars, **dset_first.analysis}
    file_vars["solution"] = config.tech.gnss_vel_comparison_report.solution.str.lower()

    # Generate figure directory to save figures generated for GNSS report
    figure_dir = config.files.path("output_gnss_vel_comparison_report_figure", file_vars=file_vars)
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    _, dfs_day, dfs_month = _generate_dataframes(dset)
    _plot_velocity_error(dfs_day, dfs_month, figure_dir, file_vars)

    # Generate GNSS comparison report
    path = config.files.path("output_gnss_vel_comparison_report", file_vars=file_vars)
    with config.files.open_path(path, create_dirs=True, mode="wt") as fid:
        rpt = Report(
            fid, rundate=dset_first.analysis["rundate"], path=path, description="Comparison of GNSS SPV analyses"
        )
        rpt.title_page()
        _add_to_report(rpt, figure_dir, dfs_day, dfs_month, file_vars)
        rpt.markdown_to_pdf()


#
# AUXILIARY FUNCTIONS
#
        
def _apply(df: pd.core.frame.DataFrame, sample: str, func: str) -> pd.core.frame.DataFrame:
    """Resample dataframe and apply given function 

    Args:
        df:      Dataframe.
        sample:  Sample definition ("D": day, "M": month)
        func:    Function to be applied ("mean", "percentile", "rms", "std")

    Returns:
        Resampled dataframe by applying given function
    """
    df_sampled = df.set_index("time_gps").resample(sample)

    if func == "mean":
        df_sampled = df_sampled.mean()
    elif func == "percentile":
        df_sampled = df_sampled.apply(lambda x: np.nanpercentile(x, q=95))
    elif func == "rms":
        df_sampled = df_sampled.apply(lambda x: np.sqrt(np.nanmean(np.square(x))))
    elif func == "std":
        df_sampled = df_sampled.std()
    else:
        log.fatal(f"Function '{func}' is not defined.")

    return df_sampled  


def _add_to_report(
    rpt: "Report",
    figure_dir: PosixPath,
    dfs_day: Dict[str, pd.core.frame.DataFrame],
    dfs_month: Dict[str, pd.core.frame.DataFrame],
    file_vars: Dict[str, Any],
) -> None:
    """Add figures and tables to report

    Args:
        rpt:         Report object.
        figure_dir:  Figure directory.
        dfs_day:     Dictionary with function type as keys ('mean', 'percentile', 'rms', 'std') and a dictionary as
                     values. The dictionary has fields as keys (e.g. site_vel_h, site_vel_3d) and the belonging 
                     dataframe as value with DAILY samples of 95th percentile and stations as columns.
        dfs_month    Dictionary with function type as keys ('mean', 'percentile', 'rms', 'std') and a dictionary as
                     values. The dictionary has fields as keys (e.g. site_vel_h, site_vel_3d) and the belonging
                     dataframe as value with MONTHLY samples of 95th percentile and stations as columns.
        file_vars:   File variables used for file and plot title naming.
    """
    
    text_def = {
            "mean": "average",
            "percentile": "95th percentile",
            "std": "standard deviation",
            "rms": "RMS",
    }

    field_def = {
            "site_vel_east":  "East site velocity component of topocentric coordinates",
            "site_vel_north":  "North site velocity component of topocentric coordinates",
            "site_vel_up":  "Up site velocity component of topocentric coordinates",
            "site_vel_h": "2D site velocity",
            "site_vel_3d": "3D site velocity",
    }

    for type_ in dfs_day.keys():
        
        for sample in config.tech.gnss_vel_comparison_report.samples.list:
    
            sample = sample.capitalize()
            rpt.add_text(f"\n# {sample} {text_def[type_]} for given solutions\n\n")
    
            if sample == "Daily":
                for field in field_def.keys():
                    dfs_day[type_][field].index = dfs_day[type_][field].index.strftime("%d-%m-%Y")
                    rpt.add_text(f"Daily {text_def[type_]} {field.upper()} results in meter/second:")
                    rpt.write_dataframe_to_markdown(dfs_day[type_][field], format="6.3f", statistic=True)
    
    
            elif sample == "Monthly":
                for field in field_def.keys():
                    rpt.add_text(f"Monthly {text_def[type_]} {field.upper()} results in meter/second:")
                    rpt.write_dataframe_to_markdown(dfs_month[type_][field], format="6.3f")
    
            # Add 2D and 3D velocity plots
            for field in field_def.keys():
                rpt.add_figure(
                    f"{figure_dir}/plot_{type_}_{field}_{sample.lower()}_{file_vars['date']}_{file_vars['solution'].lower()}.{FIGURE_FORMAT}",
                    caption=f"{text_def[type_].capitalize()} for {field_def[field]}.",
                    clearpage=True,
                )


def _generate_dataframes(dset: Dict[str, "Dataset"]) -> Dict[str, pd.core.frame.DataFrame]:
    """Generate dataframe based on station datasets

    The dataframe for each station in dictionary "dfs" has following columns:

        site_vel_h:         Horizontal site velocity
        site_vel_east:      Site velocity east component of topocentric coordinates
        site_vel_north:     Site velocity north component of topocentric coordinates
        site_vel_up:        Site velocity up component of topocentric coordinates
        site_vel_3d:        3D site velocity

    Example for "dfs" dictionary:
     
             'hons':                   time.gps     site_vel_h    site_vel_3d
                        0      2019-03-01 00:00:00    0.301738       0.057244 
                        1      2019-03-01 00:00:00    0.301738       0.057244 

             'krss':                   time.gps     site_vel_h    site_vel_3d    
                        0      2019-03-01 00:00:00    0.710014       0.186791 
                        1      2019-03-01 00:00:00    0.710014       0.186791 

    Example for "dfs_day" dictionary for "mean" key:
        'mean':{
             'site_vel_h':              nabf      vegs      hons      krss
                        time.gps                                          
                        2019-03-01  1.368875  0.935687  1.136763  0.828754
                        2019-03-02  0.924839  0.728280  0.911677  0.854832


             'site_vel_3d':             nabf      vegs      hons      krss
                        time.gps                                          
                        2019-03-01  1.715893  1.147265  1.600330  0.976541
                        2019-03-02  1.533437  1.307373  1.476295  1.136991
        }

    Example for "dfs_month" dictionary for "mean" key:
        'mean':{
            'site_vel_h':             nabf      vegs      hons      krss
                        Mar-2019  1.186240  0.861718  1.095827  1.021354
                        Apr-2019  0.891947  0.850343  0.977908  0.971099

            'site_vel_3d':            nabf      vegs      hons      krss
                        Mar-2019  1.854684  1.291406  1.450466  1.225467
                        Apr-2019  1.964404  1.706507  1.687994  1.500742
        }


    Args:
        dset: Dictionary with station name as keys and the belonging Dataset as value

    Returns:
        Tuple with following entries:

        | Element              | Description                                                                          |
        |----------------------|--------------------------------------------------------------------------------------|
        | dfs                  | Dictionary with station name as keys and the belonging dataframe as value with       |
        |                      | following dataframe columns: site_vel_h, site_vel_3d                                 |
        | dfs_day              | Dictionary with function type as keys ('mean', 'percentile', 'rms', 'std') and a     |
        |                      | dictionary as values. The dictionary has fields as keys (e.g. site_vel_h,            |
        |                      | site_vel_3d) and the belonging dataframe as value with DAILY samples of 95th         | 
        |                      | percentile and stations as columns.                                                  |
        | dfs_month            | Dictionary with function type as keys ('mean', 'percentile', 'rms', 'std') and a     |
        |                      | dictionary as values. The dictionary has fields as keys (e.g. site_vel_h,            |
        |                      | site_vel_3d) and the belonging dataframe as value with MONTHLY samples of 95th       | 
        |                      | percentile and stations as columns.                                                  |
        | dfs_month            | Dictionary with fields as keys (e.g. site_vel_h, site_vel_3d) and the belonging      |
        |                      | dataframe as value with MONTHLY samples of 95th percentile and stations as columns.  |
    """
    dsets = dset
    dfs = {}
    fields = {
        "site_vel_east": pd.DataFrame(), 
        "site_vel_north": pd.DataFrame(), 
        "site_vel_up": pd.DataFrame(), 
        "site_vel_h": pd.DataFrame(), 
        "site_vel_3d": pd.DataFrame(), 
    }
    dfs_day = { 
        "mean": copy.deepcopy(fields), 
        "percentile": copy.deepcopy(fields), 
        "std": copy.deepcopy(fields),
        "rms": copy.deepcopy(fields),
    }
    dfs_month = {
        "mean": copy.deepcopy(fields), 
        "percentile": copy.deepcopy(fields), 
        "std": copy.deepcopy(fields),
        "rms": copy.deepcopy(fields),
    }

    for station, dset in dsets.items():
        
        
        # Add necessary site velocity fields to dataset
        if "site_vel_3d" not in dset.fields:
            gnss_velocity_fields(dset)

        if dset.num_obs == 0:
            log.warn(f"Dataset '{station}' is empty.")
            continue

        # Determine dataframe with site_vel_h and site_vel_3d columns
        df = dset.as_dataframe(fields=[
                "time.gps", 
                "site_vel_east", 
                "site_vel_north", 
                "site_vel_up",
                "site_vel_h", 
                "site_vel_3d",
        ])
    
        if df.empty:
            continue
        else:
            # Save data in dictionaries
            dfs.update({station: df})

            for type_ in dfs_day.keys():
                
                df_day = _apply(df, "D", type_)
                for field in fields.keys():
                    if dfs_day[type_][field].empty:
                        dfs_day[type_][field][station] = df_day[field]
                    else:
                        dfs_day[type_][field] = pd.concat([dfs_day[type_][field], df_day[field]], axis=1)
                    dfs_day[type_][field] = dfs_day[type_][field].rename(columns={field: station})

                df_month = _apply(df, "M", type_)
                df_month.index = df_month.index.strftime("%b-%Y")
                for field in fields.keys():
                    dfs_month[type_][field][station] = df_month[field]
                    
    return dfs, dfs_day, dfs_month


def _plot_velocity_error(
    dfs_day: Dict[str, pd.core.frame.DataFrame],
    dfs_month: Dict[str, pd.core.frame.DataFrame],
    figure_dir: "pathlib.PosixPath",
    file_vars: Dict[str, Any],
) -> None:
    """Plot site velocity error plots 

    Args:
       dfs_day:     Dictionary with function type as keys ('mean', 'percentile', 'rms', 'std') and a dictionary as 
                    values. The dictionary has fields as keys (e.g. site_vel_h, site_vel_3d) and the belonging 
                    dataframe as value with DAILY samples of 95th percentile and stations as columns.  
       dfs_month:   Dictionary with function type as keys ('mean', 'percentile', 'rms', 'std') and a dictionary as 
                    values. The dictionary has fields as keys (e.g. site_vel_h, site_vel_3d) and the belonging 
                    dataframe as value with MONTHLY samples of 95th percentile and stations as columns.
       figure_dir:  Figure directory
    """
    
    ylabel_def = {
            "mean": "MEAN",
            "percentile": "95%",
            "rms": "RMS",
            "std": "STD",
    }
    
    field_def = {
        "site_vel_east":  "EAST VE",
        "site_vel_north":  "NORTH VE",
        "site_vel_up":  "UP VE",
        "site_vel_h":  "2D VE",
        "site_vel_3d":  "3D VE",
    }

    opt_args = {
        "colormap": "tab20",
        "figsize": (7, 3),
        # "grid": True,
        "marker": "o",
        "markersize": "4",
        "linestyle": "solid",
        "plot_to": "file",
        "plot_type": "plot",
        # "statistic": ["rms", "mean", "std", "min", "max", "percentile"], #TODO: Is only shown for data, which are plotted at last.
        "title": file_vars["solution"].upper(),
    }

    colors = (
        config.tech.gnss_vel_comparison_report.colors.list
        if config.tech.gnss_vel_comparison_report.colors.list
        else ["orange", "red", "violet", "blue", "green"]
    )

    colors = (
        config.tech.gnss_vel_comparison_report.colors.list
        if config.tech.gnss_vel_comparison_report.colors.list
        else ["orange", "red", "violet", "blue", "green"]
    )

    # Loop over statistical solutions
    for type_ in dfs_day.keys():

        # Get used samples
        samples = dict()
        for sample in config.tech.gnss_vel_comparison_report.samples.list:
            if "daily" == sample:
                samples["daily"] = dfs_day[type_]
            elif "monthly" == sample:
                samples["monthly"] = dfs_month[type_]
            else:
                log.fatal(f"Sample '{sample}' is not defined. Only 'daily' and/or 'monthly' can be chosen as sample.")

        # Loop over sampled data       
        for sample, sample_data in samples.items():
            
            # Loop over fields to plot
            for field in ["site_vel_east", "site_vel_north", "site_vel_up", "site_vel_h", "site_vel_3d"]:
                
                # Get y-range limits
                if field == "site_vel_h":
                    ylim = config.tech.gnss_vel_comparison_report.ylim_site_vel_h.list
                elif field == "site_vel_3d":
                    ylim = config.tech.gnss_vel_comparison_report.ylim_site_vel_3d.list
                else:
                    ylim = config.tech.gnss_vel_comparison_report.ylim.list
                    
                opt_args["ylim"] = [float(ylim[0]), float(ylim[1])] if ylim else ylim
    
                #if field == "site_vel_h":
                #    opt_args["ylim"] = [0.0, 0.03]
                #elif field == "site_vel_3d":
                #    opt_args["ylim"] = [0.0, 0.07]
    
                # Generate x- and y-arrays for plotting
                x_arrays = []
                y_arrays = []
                labels = []
                
                for station in sample_data[field].columns:
                    # if sample == "monthly":
                    #    opt_args.update({"xlim": "auto", "ylim": [0.0, 3.0]})
                    x_data = (
                        sample_data[field].index.to_pydatetime()
                        if isinstance(sample_data[field].index, pd.core.indexes.datetimes.DatetimeIndex)
                        else sample_data[field].index
                    )
                    x_arrays.append(list(x_data))
                    y_arrays.append(list(sample_data[field][station]))
                    labels.append(station.upper())
    
                # Generate plot
                plot(
                    x_arrays=x_arrays,
                    y_arrays=y_arrays,
                    xlabel="Time [GPS]",
                    ylabel=f"{field_def[field]} {ylabel_def[type_]}",
                    y_unit="m/s",
                    labels=labels,
                    colors=colors,
                    figure_path=figure_dir
                    / f"plot_{type_}_{field}_{sample}_{file_vars['date']}_{file_vars['solution'].lower()}.{FIGURE_FORMAT}",
                    opt_args=opt_args,
                )
