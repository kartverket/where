"""Compare different GNSS SPV Where datasets

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
    writers.write_one('gnss_spv_comparison_report', dset=dset, do_report=False)

"""
# Standard library imports
from typing import Any, Dict

# External library imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Midgard imports
from midgard.dev import plugins
from midgard.plot.matplotlib_extension import plot

# Where imports
import where
from where.lib import config
from where.lib import log
from where.writers._report import Report

FIGURE_FORMAT = "png"
FILE_NAME = __name__.split(".")[-1]


@plugins.register
def gnss_spv_comparison_report(dset: Dict[str, "Dataset"]) -> None:
    """Compare GNSS SPV datasets

    Args:
        dset:  Dictionary with station name as keys and the belonging Dataset as value
    """
    dset_first = dset[list(dset.keys())[0]]
    dset_first.vars["solution"] = config.tech.gnss_spv_comparison_report.solution.str.lower()

    # Generate figure directory to save figures generated for GNSS report
    figure_dir = config.files.path("output_gnss_spv_comparison_report_figure", file_vars=dset_first.vars)
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    _, dfs_day, dfs_month = _generate_dataframes(dset)
    _plot_velocity_error(dfs_day, dfs_month, figure_dir, dset_first.vars)

    # Generate GNSS comparison report
    path = config.files.path("output_gnss_spv_comparison_report", file_vars=dset_first.vars)
    with config.files.open_path(path, create_dirs=True, mode="wt") as fid:
        rpt = Report(fid, rundate=dset_first.rundate, path=path, description="Comparison of GNSS SPV analyses")
        rpt.title_page()
        _add_to_report(rpt, figure_dir, dfs_day, dfs_month, dset_first.vars)
        rpt.markdown_to_pdf()


def _add_to_report(
    rpt: "Report",
    figure_dir: "pathlib.PosixPath",
    dfs_day: Dict[str, pd.core.frame.DataFrame],
    dfs_month: Dict[str, pd.core.frame.DataFrame],
    file_vars: Dict[str, Any],
) -> None:
    """Add figures and tables to report

    Args:
        rpt:         Report object.
        figure_dir:  Figure directory.
        dfs_day:     Dictionary with fields as keys (e.g. 2d_vel, 3d_vel) and the belonging dataframe as value with DAILY
                     samples of 95th percentile and stations as columns.
        dfs_month:   Dictionary with fields as keys (e.g. 2d_vel, 3d_vel) and the belonging dataframe as value with MONTHLY
                     samples of 95th percentile and stations as columns.
        file_vars:   File variables used for file and plot title naming.
    """

    for sample_name in ["Daily", "Monthly"]:

        rpt.add_text(f"\n# {sample_name} 95th percentile 2D and 3D solutions\n\n")

        if sample_name == "Daily":
            for field in dfs_day.keys():
                dfs_day[field].index = dfs_day[field].index.strftime("%d-%m-%Y")

            rpt.add_text("Daily 95th percentile 2D velocity results in meter/second:")
            rpt.write_dataframe_to_markdown(dfs_day["hpe"], format="6.2f", statistic=True)

            rpt.add_text("Daily 95th percentile 3D velocity in meter/second:")
            rpt.write_dataframe_to_markdown(dfs_day["vpe"], format="6.2f", statistic=True)

        elif sample_name == "Monthly":
            rpt.add_text("Monthly 95th percentile 2D velocity results in meter/second:")
            rpt.write_dataframe_to_markdown(dfs_month["hpe"], format="6.2f")

            rpt.add_text("Monthly 95th percentile 3D velocity results in meter/second:")
            rpt.write_dataframe_to_markdown(dfs_month["vpe"], format="6.2f")

        # Add 2D and 3D velocity plots
        rpt.add_figure(
            f"{figure_dir}/plot_hpe_{sample_name.lower()}_{file_vars['date']}_{file_vars['solution'].lower()}.{FIGURE_FORMAT}",
            caption="95th percentile for 2D velocity.",
            clearpage=True,
        )

        rpt.add_figure(
            f"{figure_dir}/plot_vpe_{sample_name.lower()}_{file_vars['date']}_{file_vars['solution'].lower()}.{FIGURE_FORMAT}",
            caption="95th percentile for 3D velocity.",
            clearpage=True,
        )


def _generate_dataframes(dset: Dict[str, "Dataset"]) -> Dict[str, pd.core.frame.DataFrame]:
    """Generate dataframe based on station datasets

    The dataframe for each station in dictionary "dfs" has following columns:

        east:  East-coordinate in topocentric system
        north: North-coordinate in topocentric system
        up:    Up-coordinate in topocentric system
        hpe:   horizontal position error
        vpe:   vertical position error

    Example for "dfs" dictionary:
     
             'hons':                   time.gps       2d_vel    3d_vel 
                        0      2019-03-01 00:00:00  0.301738  0.057244 
                        1      2019-03-01 00:00:00  0.301738  0.057244 

             'krss':                   time.gps       2d_vel    3d_vel    
                        0      2019-03-01 00:00:00  0.710014  0.186791 
                        1      2019-03-01 00:00:00  0.710014  0.186791 

    Example for "dfs_day" dictionary:

             '2d_vel':              nabf      vegs      hons      krss
                        time.gps                                          
                        2019-03-01  1.368875  0.935687  1.136763  0.828754
                        2019-03-02  0.924839  0.728280  0.911677  0.854832


             '3d_vel':              nabf      vegs      hons      krss
                        time.gps                                          
                        2019-03-01  1.715893  1.147265  1.600330  0.976541
                        2019-03-02  1.533437  1.307373  1.476295  1.136991

    Example for "dfs_month" dictionary:

            '2d_vel':             nabf      vegs      hons      krss
                        Mar-2019  1.186240  0.861718  1.095827  1.021354
                        Apr-2019  0.891947  0.850343  0.977908  0.971099

            '3d_vel':             nabf      vegs      hons      krss
                        Mar-2019  1.854684  1.291406  1.450466  1.225467
                        Apr-2019  1.964404  1.706507  1.687994  1.500742


    Args:
        dset: Dictionary with station name as keys and the belonging Dataset as value

    Returns:
        Tuple with following entries:

        | Element              | Description                                                                          |
        |----------------------|--------------------------------------------------------------------------------------|
        | dfs                  | Dictionary with station name as keys and the belonging dataframe as value with       |
        |                      | following dataframe columns: 2d_vel, 3d_vel                                          |
        | dfs_day              | Dictionary with fields as keys (e.g. 2d_vel, 3d_vel) and the belonging dataframe as  |
        |                      | value with DAILY samples of 95th percentile and stations as columns.                 |
        | dfs_month            | Dictionary with fields as keys (e.g. 2d_vel, 3d_vel) and the belonging dataframe as  |
        |                      | value with MONTHLY samples of 95th percentile and stations as columns.               |
    """
    dsets = dset
    dfs = {}
    dfs_day = {"2d_vel": pd.DataFrame(), "3d_vel": pd.DataFrame()}
    dfs_month = {"2d_vel": pd.DataFrame(), "3d_vel": pd.DataFrame()}

    for station, dset in dsets.items():

        if dset.num_obs == 0:
            log.warn(f"Dataset '{station}' is empty.")
            continue

        # Determine dataframe with 2d_vel and 3d_vel columns
        # TODO: How to ensure that GPS time scale is used? fields=["time.gps", ...] does not work longer.
        df = dset.as_dataframe(fields=["time", "2d_vel", "3d_vel"])
        if df.empty:
            continue
        else:
            # Save data in dictionaries
            dfs.update({station: df})

            # TODO This does not work anymore for Pandas version 1.0: df_day = df.set_index("time").resample("D", how=lambda x: np.nanpercentile(x, q=95))
            df_day = df.set_index("time").resample("D").apply(lambda x: np.nanpercentile(x, q=95))
            for field in dfs_day.keys():
                if dfs_day[field].empty:
                    dfs_day[field][station] = df_day[field]
                else:
                    dfs_day[field] = pd.concat([dfs_day[field], df_day[field]], axis=1)
                dfs_day[field] = dfs_day[field].rename(columns={field: station})

            # TODO This does not work anymore for Pandas version 1.0: df_month = df.set_index("time").resample("M", how=lambda x: np.nanpercentile(x, q=95))
            df_month = df.set_index("time").resample("M").apply(lambda x: np.nanpercentile(x, q=95))
            df_month.index = df_month.index.strftime("%b-%Y")
            for field in dfs_month.keys():
                dfs_month[field][station] = df_month[field]

    return dfs, dfs_day, dfs_month


#
# PLOT FUNCTIONS
#
def _plot_velocity_error(
    dfs_day: Dict[str, pd.core.frame.DataFrame],
    dfs_month: Dict[str, pd.core.frame.DataFrame],
    figure_dir: "pathlib.PosixPath",
    file_vars: Dict[str, Any],
) -> None:
    """Plot 2D and 3D velocity error plots (95th percentile) 

    Args:
       dfs_day:     Dictionary with fields as keys (e.g. 2d_vel, 3d_vel) and the belonging dataframe as value with DAILY
                    samples of 95th percentile and stations as columns.
       dfs_month:   Dictionary with fields as keys (e.g. hpe, vpe) and the belonging dataframe as value with MONTHLY
                    samples of 95th percentile and stations as columns.
       figure_dir:  Figure directory
    """
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
        config.tech.gnss_spv_comparison_report.colors.list
        if config.tech.gnss_spv_comparison_report.colors.list
        else ["orange", "red", "violet", "blue", "green"]
    )

    colors = (
        config.tech.gnss_spv_comparison_report.colors.list
        if config.tech.gnss_spv_comparison_report.colors.list
        else ["orange", "red", "violet", "blue", "green"]
    )

    # Loop over sampled data
    samples = {"daily": dfs_day, "monthly": dfs_month}
    for sample_name, sample_data in samples.items():

        # Loop over fields to plot
        for field in ["2d_vel", "3d_vel"]:

            if field == "2d_vel":
                opt_args["ylim"] = [0.0, 1.0]
            elif field == "3d_vel":
                opt_args["ylim"] = [0.0, 1.0]

            # Generate x- and y-arrays for plotting
            x_arrays = []
            y_arrays = []
            labels = []

            for station in sample_data[field].columns:
                if sample_name == "monthly":
                    opt_args.update({"xlim": "auto", "ylim": [0.0, 3.0]})
                x_arrays.append(list(sample_data[field].index))
                y_arrays.append(list(sample_data[field][station]))
                labels.append(station.upper())

            # Generate plot
            plot(
                x_arrays=x_arrays,
                y_arrays=y_arrays,
                xlabel="Time [GPS]",
                ylabel=f"{field.upper()} 95%",
                y_unit="m",
                labels=labels,
                colors=colors,
                figure_path=figure_dir
                / f"plot_{field}_{sample_name}_{file_vars['date']}_{file_vars['solution'].lower()}.{FIGURE_FORMAT}",
                opt_args=opt_args,
            )
