"""Write report about a GNSS analysis run

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
from midgard.plot.matplotlib_extension import plot_scatter_subplots, plot

# Where imports
from where.data import dataset3 as dataset
from where.lib import config
from where.lib import log
from where.data import position
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
        "gnss_range", 
        "delay", 
        "Correction of range between satellite and receiver",
    ),
    PlotField(
        "gnss_satellite_clock",
        "delay",
        "Correction of satellite clock",
    ),
    PlotField(
        "troposphere_radio",
        "delay",
        "Correction of tropospheric delay",
    ),
    PlotField(
        "gnss_total_group_delay",
        "delay",
        "Correction of total group delay",
    ),
    PlotField(
        "gnss_ionosphere",
        "delay",
        "Correction of ionospheric delay",
    ),
    PlotField(
        "gnss_relativistic_clock",
        "delay",
        "Correction of relativistic clock effect due to orbit eccentricity",
    ),
   # PlotField(
   #     "estimate_gnss_rcv_clock",
   #     ("estimate_gnss_rcv_clock",),
   #     "Estimate of receiver clock",
   # ),
)


@plugins.register
def gnss_report(dset: "Dataset") -> None:
    """Write report about a GNSS analysis run

    Args:
        dset:        A dataset containing the data.
    """
    file_vars={**dset.vars, **dset.analysis}

    # TODO: Better solution?
    if "station" not in file_vars:  # necessary if called for example by ./where/tools/concatenate.py
        file_vars["station"] = ""
        file_vars["STATION"] = ""

    # Generate figure directory to save figures generated for GNSS report
    figure_dir = config.files.path("output_gnss_report_figure", file_vars=file_vars)
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    if config.tech.gnss_report.kinematic.bool:
        _plot_position_kinematic(dset, figure_dir)
    else:
        _plot_position(dset, figure_dir)
    _plot_residual(dset, figure_dir)

    # Generate GNSS report
    path = config.files.path("output_gnss_report", file_vars=file_vars)
    with config.files.open_path(path, create_dirs=True, mode="wt") as fid:
        rpt = Report(fid, rundate=dset.analysis["rundate"], path=path, description="GNSS analysis")
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
    rpt.add_text("\n# GNSS site position analysis\n\n")

    if config.tech.gnss_report.kinematic.bool:

        # Plot site position 
        rpt.add_figure(
            f"{figure_dir}/plot_timeseries_llh.{FIGURE_FORMAT}",
            caption="Site position given in geodetic coordinates (latitude, longitude, height).",
            clearpage=True,
        )

        # Plot horizontal error
        rpt.add_figure(
            f"{figure_dir}/plot_horizontal_position.{FIGURE_FORMAT}",
            caption="Horizontal position",
            clearpage=True,
        )

    else:
        # Plot site position vs. reference position timeseries
        rpt.add_figure(
            f"{figure_dir}/plot_timeseries_enu.{FIGURE_FORMAT}",
            caption="Site position vs. reference position given in topocentric coordinates (East, North, Up).",
            clearpage=True,
        )

        # Plot horizontal error
        rpt.add_figure(
            f"{figure_dir}/plot_horizontal_error.{FIGURE_FORMAT}",
            caption="Horizontal error of site position vs. reference position",
            clearpage=True,
        )

        # Plot HPE and VPE timeseries
        rpt.add_figure(
            f"{figure_dir}/plot_timeseries_pdop_hpe_vpe.{FIGURE_FORMAT}",
            caption="Horizontal and vertical position error of site position in relation to reference position.",
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
                caption=f"Skyplot for {enums.gnss_id_to_name[figure_path.stem[-1]]}", 
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
def _plot_position(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Plot site position plots

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """

    # TODO: Should be replaced by dset.site_pos_vs_ref_east, ... (see gnss_position.py writer)
    ref_pos = position.Position(
        val=np.array([dset.meta["pos_x"], dset.meta["pos_y"], dset.meta["pos_z"]]), system="trs"
    )
    enu = (dset.site_pos.trs.pos - ref_pos).enu

    figure_path = figure_dir / f"plot_timeseries_enu.{FIGURE_FORMAT}"
    log.debug(f"Plot {figure_path}.")
    plot_scatter_subplots(
        x_array=dset.time.gps.datetime,
        y_arrays=[enu.east, enu.north, enu.up],
        xlabel="Time [GPS]",
        ylabels=["East", "North", "Up"],
        colors=["steelblue", "darkorange", "limegreen"],
        y_units=["meter", "meter", "meter"],
        figure_path=figure_path,
        opt_args={
            "figsize": (6, 6.8),
            "plot_to": "file",
            "sharey": True,
            "title": "Site position vs. reference position",
            "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
        },
    )

    # TODO: Should be replaced by dset.hpe, dset.vpe (see gnss_position.py writer)
    hpe = np.sqrt(enu.east ** 2 + enu.north ** 2)
    vpe = np.absolute(enu.up)

    figure_path = figure_dir / f"plot_timeseries_pdop_hpe_vpe.{FIGURE_FORMAT}"
    log.debug(f"Plot {figure_path}.")
    plot_scatter_subplots(
        x_array=dset.time.gps.datetime,
        y_arrays=[dset.pdop, hpe, vpe],
        xlabel="Time [GPS]",
        ylabels=["PDOP", "HPE", "VPE"],
        colors=["steelblue", "darkorange", "limegreen"],
        y_units=[None, "meter", "meter"],
        figure_path=figure_path,
        opt_args={
            "figsize": (6, 5),
            "plot_to": "file",
            "sharey": False,
            # "title": "Horizontal and vertical position error vs. reference position",
            "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
        },
    )

    figure_path = figure_dir / f"plot_horizontal_error.{FIGURE_FORMAT}"
    log.debug(f"Plot {figure_path}.")
    plot_scatter_subplots(
        x_array=enu.east,
        y_arrays=[enu.north],
        xlabel="East [meter]",
        ylabels=["North"],
        y_units=["meter"],
        figure_path=figure_path,
        opt_args={
            "grid": True,
            "figsize": (6, 6),
            "histogram": "x, y",
            "plot_to": "file",
            "title": "Horizontal error",
            "xlim": [-8, 8],
            "ylim": [-8, 8],
        },
    )


def _plot_position_kinematic(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Plot site position plots for kinematic solution

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """

    figure_path = figure_dir / f"plot_timeseries_llh.{FIGURE_FORMAT}"
    log.debug(f"Plot {figure_path}.")
    plot_scatter_subplots(
        x_array=dset.time.gps.datetime,
        y_arrays=[np.rad2deg(dset.site_pos.llh.lat), np.rad2deg(dset.site_pos.llh.lon), dset.site_pos.llh.height],
        xlabel="Time [GPS]",
        ylabels=["Latitude", "Longitude", "Height"],
        colors=["steelblue", "darkorange", "limegreen"],
        y_units=["degree", "degree", "meter"],
        figure_path=figure_path,
        opt_args={
            "figsize": (6, 6.8),
            "plot_to": "file",
            "sharey": False,
            "title": "Site position",
            "statistic": ["min", "max"],
        },
    )


    figure_path = figure_dir / f"plot_horizontal_position.{FIGURE_FORMAT}"
    log.debug(f"Plot {figure_path}.")
    lon = np.rad2deg(dset.site_pos.llh.lon)
    lat = np.rad2deg(dset.site_pos.llh.lat)

    # Determine range
    dx = np.abs(np.max(lon) - np.min(lon))
    dy = np.abs(np.max(lat) - np.min(lat))
    incr = np.max([dx, dy])/2 # increment

    plot_scatter_subplots(
        x_array=lon,
        y_arrays=[lat],
        xlabel="Longitude [degree]",
        ylabels=["Latitude"],
        y_units=["degree"],
        figure_path=figure_path,
        opt_args={
            "grid": True,
            "figsize": (6, 6),
            "plot_to": "file",
            "title": "Horizontal position",
            "xlim": [np.mean(lon) - incr, np.mean(lon) + incr],
            "ylim": [np.mean(lat) - incr, np.mean(lat) + incr],
        },
    )

    # TODO: Determine speed - this is wrong
    #for epoch in dset.time.gps.gps_seconds:
    #    idx = dset.time.gps.gps_seconds == epoch
    #
    #    dx = np.insert(np.diff(dset.site_pos.trs.x[idx]), 0, float('nan'))
    #    dy = np.insert(np.diff(dset.site_pos.trs.y[idx]), 0, float('nan'))
    #    dz = np.insert(np.diff(dset.site_pos.trs.z[idx]), 0, float('nan'))
    #    dt = np.insert(np.diff(dset.time.gps.gps_seconds[idx]), 0, float('nan'))
    #    speed = np.sqrt(np.power(dx, 2) * np.power(dy, 2) * np.power(dz, 2)) / dt
        


def _plot_residual(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Plot residual plot

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """
    figure_path = figure_dir / f"plot_residual.{FIGURE_FORMAT}"
    log.debug(f"Plot {figure_path}")

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

    plot(
        x_arrays=x_arrays,
        y_arrays=y_arrays,
        xlabel="Time [GPS]",
        ylabel="Post-fit residual",
        y_unit="meter",
        colors=colors,
        figure_path=figure_path,
        opt_args={
            "figsize": (7, 4),
            "histogram": "y",
            "histogram_size": 0.8,
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

