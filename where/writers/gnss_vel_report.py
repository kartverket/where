"""Write report about a GNSS site velocity analysis run

Description:
------------



"""
# Standard library imports
from enum import Enum
from collections import namedtuple
from typing import Union

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.collections import enums 
from midgard.dev import plugins
from midgard.gnss import gnss
from midgard.plot.matplotext import MatPlotExt
from midgard.math import rotation

# Where imports
from where.data import dataset3 as dataset
from where.lib import config
from where.lib import log
from where.writers._gnss_plot import GnssPlot
from where.writers._report import Report


FIGURE_DPI = 200
FIGURE_FORMAT = "png"

PlotField = namedtuple("PlotField", ["name", "collection", "caption"])
PlotField.__new__.__defaults__ = (None,) * len(PlotField._fields)
PlotField.__doc__ = """A convenience class for defining a output field for plotting

    Args:
        name  (str):              Unique name
        collection (Tuple[str]):  Collection name
        caption (str):            Caption of plot
    """

FIELDS = (
    PlotField(
        "gnss_range_rate",
        "delay",
        "Correction of range between satellite and receiver",
    ),
    PlotField(
        "gnss_satellite_clock_rate",
        "delay",
        "Correction of satellite clock rate",
    ),
    PlotField(
        "gnss_earth_rotation_drift",
        "delay",
        "Correction of Earth rotation drift",
    ),
    PlotField(
        "gnss_relativistic_clock_rate",
        "delay",
        "Correction of relativistic clock rate effect due to orbit eccentricity",
    ),
   # PlotField(
   #     "estimate_gnss_rcv_clock_rate",
   #     "estimate_gnss_rcv_clock_rate",
   #     "Estimate of receiver clock rate",
   # ),
)


@plugins.register
def gnss_vel_report(dset: "Dataset") -> None:
    """Write report about a GNSS velocity analysis run

    Args:
        dset:        A dataset containing the data.
    """
    file_vars = {**dset.vars, **dset.analysis}
    # TODO: Better solution?
    if "station" not in file_vars:  # necessary if called for example by ./where/tools/concatenate.py
        file_vars["station"] = ""
        file_vars["STATION"] = ""

    # Generate figure directory to save figures generated for GNSS report
    figure_dir = config.files.path("output_gnss_vel_report_figure", file_vars=file_vars)
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    _plot_velocity(dset, figure_dir)
    _plot_residual(dset, figure_dir)

    # Generate GNSS velocity report
    path = config.files.path("output_gnss_vel_report", file_vars=file_vars)
    with config.files.open_path(path, create_dirs=True, mode="wt") as fid:
        rpt = Report(fid, rundate=dset.analysis["rundate"], path=path, description="GNSS site velocity analysis")
        rpt.title_page()
        rpt.write_config()
        _add_to_report(dset, rpt, figure_dir)
        rpt.markdown_to_pdf()


def _add_to_report(dset: "Dataset", rpt: "Report", figure_dir: "pathlib.PosixPath") -> None:
    """Add figures and tables to report

    Args:
        dset:        A dataset containing the data.
        rpt:         Report object.
        figure_dir:  Figure directory.
    """
    plt = GnssPlot(dset, figure_dir)

    #
    # Position
    #
    rpt.add_text("\n# GNSS site velocity analysis\n\n")

    # Plot site velocity
    rpt.add_figure(
        f"{figure_dir}/plot_timeseries_enu.{FIGURE_FORMAT}",
        caption="Site velocity in topocentric coordinates (East, North, Up).",
        clearpage=True,
    )

    # Plot horizontal error
    rpt.add_figure(
        f"{figure_dir}/plot_horizontal_velocity.{FIGURE_FORMAT}",
        caption="Horizontal velocity",
        clearpage=True,
    )

    # Plot 3D timeseries
    rpt.add_figure(
        f"{figure_dir}/plot_timeseries_pdop_hv_3d.{FIGURE_FORMAT}",
        caption="PDOP together with Horizontal, vertical and 3D velocity of site position",
        clearpage=True,
    )

    #
    # Residual
    #
    rpt.add_text("\n# GNSS residual\n\n")

    # Add outlier table
    # MURKS: does not work at the moment. complement_with is not implemented in Dataset v3.
    # MURKS rpt.write_dataframe_to_markdown(_table_outlier_overview(dset))

    # Plot residuals
    rpt.add_figure(
        f"{figure_dir}/plot_residual.{FIGURE_FORMAT}",
        # MURKScaption="Post-fit residuals, whereby the red dots represents the rejected outliers. The histogram represent only number of residuals from kept observations.",
        caption="Post-fit residuals.",
        clearpage=True,
    )

    #
    # Dilution of precision (DOP)
    #
    if "pdop" in dset.fields:
        rpt.add_text("\n# Dilution of precision\n\n")

        # Plot DOP
        rpt.add_figure(
                figure_path=plt.plot_dop(), 
                caption="Dilution of precision.", 
                clearpage=True,
        )

    #
    # Satellite plots
    #
    rpt.add_text("\n# Satellite plots\n\n")

    rpt.add_figure(
        figure_path=plt.plot_number_of_satellites_used(),
        caption="Number of satellites for each observation epoch",
        clearpage=False,
    )

    figure_path = plt.plot_satellite_overview()
    if figure_path is not None:  # Note: Does not exists for concatenated Datasets.
        rpt.add_figure(
            figure_path=figure_path,
            caption="Overview over satellite observations. Red coloured: Observation rejected in orbit stage (e.g. unhealthy satellites, exceeding validity length, no orbit data available); Orange coloured: Observation rejected in edit stage; Green coloured: Kept observations after edit stage.",
            clearpage=False,
        )

    for figure_path in plt.plot_skyplot():
        rpt.add_figure(
            figure_path=figure_path,
            caption="Skyplot for {enums.gnss_id_to_name[figure_path.stem[-1]]}",
            clearpage=False,
        )

    for figure_path in plt.plot_satellite_elevation():
        rpt.add_figure(
                figure_path=figure_path,
                caption=f"Satellite elevation for {enums.gnss_id_to_name[figure_path.stem[-1]]}", 
                clearpage=True,
        )

    #
    # Model parameter plots
    #
    rpt.add_text("\n# Plots of model parameters\n\n")

    for field in FIELDS:
        
        if f"{field.collection}.{field.name}" in dset.fields:       
            for figure_path in plt.plot_field(field.name, field.collection):
                system, _ = figure_path.stem.split("_")[2:4] 
                rpt.add_figure(
                    figure_path=figure_path, 
                    caption=f"{field.caption} for {enums.gnss_id_to_name[system].value} observation", 
                    clearpage=True,
                )


#
# PLOT FUNCTIONS
#
def _plot_velocity(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Plot site velocity plots

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """

    
    lat, lon, height = dset.site_pos.pos.llh.T
    vel_enu = np.squeeze(rotation.trs2enu(lat, lon) @  dset.site_vel[:,:,None]) 

    plt = MatPlotExt()
    plt.plot_subplots(
        x_array=dset.time.gps.datetime,
        y_arrays=[vel_enu[:, 0], vel_enu[:, 1], vel_enu[:, 2]],
        xlabel="Time [GPS]",
        ylabels=["East", "North", "Up"],
        colors=["steelblue", "darkorange", "limegreen"],
        y_units=["m/s", "m/s", "m/s"],
        figure_path=figure_dir / f"plot_timeseries_enu.{FIGURE_FORMAT}",
        options={
            "figsize": (6, 6.8),
            "plot_to": "file",
            "sharey": True,
            "title": "Site velocity",
            "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
        },
    )

    vel_h = np.sqrt(vel_enu[:,0] ** 2 + vel_enu[:,1] ** 2) 
    vel_v = np.absolute(vel_enu[:,2])
    #vel_3d = np.sqrt(vel_enu[:,0] ** 2 + vel_enu[:,1] ** 2 + vel_enu[:,2] ** 2)
    vel_3d = np.sqrt(dset.site_vel[:,0] ** 2 + dset.site_vel[:,1] ** 2 + dset.site_vel[:,2] ** 2)

    plt = MatPlotExt()
    plt.plot_subplots(
        x_array=dset.time.gps.datetime,
        y_arrays=[dset.pdop, vel_h, vel_v, vel_3d],
        xlabel="Time [GPS]",
        ylabels=["PDOP", "HV", "VV", "3D"],
        colors=["steelblue", "darkorange", "limegreen", "red"],
        y_units=[None, "m/s", "m/s", "m/s"],
        figure_path=figure_dir / f"plot_timeseries_pdop_hv_3d.{FIGURE_FORMAT}",
        options={
            "figsize": (7, 7),
            "plot_to": "file",
            "sharey": False,
            # "title": "2D (horizontal) and 3D velocity",
            "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
        },
    )

    plt = MatPlotExt()
    plt.plot_subplots(
        x_array=vel_enu[:, 0],
        y_arrays=[vel_enu[:, 1]],
        xlabel="East [m/s]",
        ylabels=["North"],
        y_units=["m/s"],
        figure_path=figure_dir / f"plot_horizontal_velocity.{FIGURE_FORMAT}",
        options={
            "grid": True,
            "figsize": (6, 6),
            "histogram": "x, y",
            "histogram_binwidth": 0.002,
            "plot_to": "file",
            "title": "Horizontal velocity",
            "xlim": [-0.1, 0.1],
            "ylim": [-0.1, 0.1],
        },
    )


def _plot_residual(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Plot residual plot

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """
    figure_path = figure_dir / f"plot_residual.{FIGURE_FORMAT}"
    dset_outlier = _get_outliers_dataset(dset)

    if dset_outlier == enums.ExitStatus.error:
        # NOTE: This is the case for concatencated Datasets, where "calculate" stage data are not available.
        log.warn(f"No data for calculate stage available. No outliers are plotted in {figure_path}.")
        x_arrays = [dset.time.gps.datetime]
        y_arrays = [dset.residual]
        colors = ["dodgerblue"]
    else:
        if dset_outlier.num_obs:
            x_arrays = [dset_outlier.time.gps.datetime, dset.time.gps.datetime]
            y_arrays = [dset_outlier.residual, dset.residual]
            colors = ["red", "dodgerblue"]
        else:
            log.debug("No outliers detected.")
            x_arrays = [dset.time.gps.datetime]
            y_arrays = [dset.residual]
            colors = ["dodgerblue"]

    plt = MatPlotExt()
    plt.plot(
        x_arrays=x_arrays,
        y_arrays=y_arrays,
        xlabel="Time [GPS]",
        ylabel="Post-fit residual",
        y_unit="m/s",
        colors=colors,
        figure_path=figure_path,
        options={
            "figsize": (7, 4),
            "histogram": "y",
            "histogram_size": 0.8,
            "histogram_binwidth": 0.002,
            "plot_to": "file",
            "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
        },
    )


#
# TABLE GENERATION FUNCTIONS
#
def _table_outlier_overview(dset: "Dataset"):
    """Generate Dataframe table with overview over number of navigation messages

    Args:
        dset:      A dataset containing the data.

    Returns:
        Dataframe with satellites as indices and following columns:

        | Name        | Description                                                                                  |
        |-------------|----------------------------------------------------------------------------------------------|
        | outlier     | Number of outliers for each satellite                                                        |


        Example:

            |    |outlier | 
            |----|--------|
            | G01|      0 |
            | G02|     11 |
            | G03|      3 |
            | .. |    ... |
            | SUM|     42 |

    """
    columns = ["outlier"]
    df = pd.DataFrame(columns=columns)

    dset_outlier = _get_outliers_dataset(dset)
    if dset_outlier == enums.ExitStatus.error:
        # NOTE: This is the case for concatencated Datasets, where "calculate" stage data are not available.
        log.warn(f"No data for calculate stage available. Outliers can not be detected.")
        return df

    if dset_outlier.num_obs:
        log.debug("No outlier detected.")
        return df

    for satellite in sorted(dset.unique("satellite")):
        idx = dset_outlier.filter(satellite=satellite)
        row = [len(dset_outlier.satellite[idx])]
        df = df.append(pd.DataFrame([row], columns=columns, index=[satellite]))

    df = df.append(pd.DataFrame([[len(dset_outlier.satellite)]], columns=columns, index=["**SUM**"]))

    return df


#
# AUXILIARY FUNCTIONS
#
def _get_outliers_dataset(dset: "Dataset") -> Union["Dataset", Enum]:
    """Get dataset with outliers

    Args:
       dset:        A dataset containing the data.

    Returns:
       Dataset with outliers or error exit status if no data for "calculate" stage are available
    """

    # Get Dataset where no outliers are rejected
    file_vars = {**dset.vars, **dset.analysis}
    file_vars["stage"] = "calculate"

    try:
        dset_complete = dataset.Dataset.read(**file_vars)
    except OSError:
        log.warn(f"Could not read dataset {config.files.path('dataset', file_vars=file_vars)}.")
        return enums.ExitStatus.error

    # Get relative complement, which corresponds to "outlier" dataset
    # dset_outliers = dset_complete.complement_with(dset, complement_by=["time", "satellite"])
    dset_outliers = dataset.Dataset(num_obs=0)  # MURKS: complement_with does not exists so far in Dataset v3.

    return dset_outliers

