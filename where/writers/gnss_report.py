"""Write report about a GNSS analysis run

Description:
------------



"""
# Standard library imports
from datetime import datetime
from typing import Any, Dict, List, Tuple, Union

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from where.data import position
from midgard.dev import plugins
from midgard.plot.matplotlib_extension import plot_scatter_subplots, plot

# Where imports
from where.data import dataset3 as dataset
from where.lib import config
from where.lib import gnss
from where.lib import log
from where.writers._report import Report


FIGURE_DPI = 200
FIGURE_FORMAT = "png"


@plugins.register
def gnss_report(dset: "Dataset") -> None:
    """Write report about a GNSS analysis run

    Args:
        dset:        A dataset containing the data.
    """

    # TODO: Better solution?
    if "station" not in dset.vars:  # necessary if called for example by ./where/tools/concatenate.py
        dset.vars["station"] = ""
        dset.vars["STATION"] = ""

    # Generate figure directory to save figures generated for GNSS report
    figure_dir = config.files.path("output_gnss_report_figure", file_vars=dset.vars)
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    _plot_position(dset, figure_dir)
    _plot_residual(dset, figure_dir)
    _plot_number_of_satellites(dset, figure_dir)
    _plot_satellite_overview(dset, figure_dir)
    _plot_skyplot(dset, figure_dir)
    _plot_satellite_elevation(dset, figure_dir)

    if "pdop" in dset.fields:
        _plot_dop(dset, figure_dir)

    # Generate GNSS report
    path = config.files.path("output_gnss_report", file_vars=dset.vars)
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

    #
    # Position
    #
    rpt.add_text("\n# GNSS site position analysis\n\n")

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
    rpt.write_dataframe_to_markdown(_table_outlier_overview(dset))

    # Plot residuals
    rpt.add_figure(
        f"{figure_dir}/plot_residual.{FIGURE_FORMAT}",
        caption="Post-fit residuals, whereby the red dots represents the rejected outliers. The histogram represent only number of residuals from kept observations.",
        clearpage=True,
    )

    #
    # Dilution of precision (DOP)
    #
    if "pdop" in dset.fields:
        rpt.add_text("\n# Dilution of precision\n\n")

        # Plot DOP
        rpt.add_figure(f"{figure_dir}/plot_dop.{FIGURE_FORMAT}", caption="Dilution of precision.", clearpage=True)

    #
    # Satellite plots
    #
    rpt.add_text("\n# Satellite plots\n\n")

    rpt.add_figure(
        f"{figure_dir}/plot_number_of_satellites.{FIGURE_FORMAT}",
        caption="Number of satellites for each observation epoch",
        clearpage=False,
    )

    figure_path = figure_path = figure_dir / f"plot_satellite_overview.{FIGURE_FORMAT}"
    if figure_path.exists():  # Note: Does not exists for concatenated Datasets.
        rpt.add_figure(
            figure_path,
            caption="Overview over satellite observations. Red coloured: Observation rejected in orbit stage (e.g. unhealthy satellites, exceeding validity length, no orbit data available); Orange coloured: Observation rejected in edit stage; Green coloured: Kept observations after edit stage.",
            clearpage=False,
        )

    rpt.add_figure(f"{figure_dir}/plot_skyplot.{FIGURE_FORMAT}", caption="Skyplot", clearpage=False)

    rpt.add_figure(
        f"{figure_dir}/plot_satellite_elevation.{FIGURE_FORMAT}", caption="Satellite elevation", clearpage=True
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

    plot_scatter_subplots(
        x_array=dset.time.gps.datetime,
        y_arrays=[enu.east, enu.north, enu.up],
        xlabel="Time [GPS]",
        ylabels=["East", "North", "Up"],
        colors=["steelblue", "darkorange", "limegreen"],
        y_units=["meter", "meter", "meter"],
        figure_path=figure_dir / f"plot_timeseries_enu.{FIGURE_FORMAT}",
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

    plot_scatter_subplots(
        x_array=dset.time.gps.datetime,
        y_arrays=[dset.pdop, hpe, vpe],
        xlabel="Time [GPS]",
        ylabels=["PDOP", "HPE", "VPE"],
        colors=["steelblue", "darkorange", "limegreen"],
        y_units=[None, "meter", "meter"],
        figure_path=figure_dir / f"plot_timeseries_pdop_hpe_vpe.{FIGURE_FORMAT}",
        opt_args={
            "figsize": (6, 5),
            "plot_to": "file",
            "sharey": False,
            # "title": "Horizontal and vertical position error vs. reference position",
            "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
        },
    )

    plot_scatter_subplots(
        x_array=enu.east,
        y_arrays=[enu.north],
        xlabel="East [meter]",
        ylabels=["North"],
        y_units=["meter"],
        figure_path=figure_dir / f"plot_horizontal_error.{FIGURE_FORMAT}",
        opt_args={
            "grid": True,
            "figsize": (6, 6),
            "histogram": "x, y",
            "plot_to": "file",
            "title": "Horizontal error",
            "xlim": [-4, 4],
            "ylim": [-4, 4],
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

    if dset_outlier == 1:
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


def _plot_dop(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Plot DOP

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """
    plot(
        x_arrays=[
            dset.time.gps.datetime,
            dset.time.gps.datetime,
            dset.time.gps.datetime,
            dset.time.gps.datetime,
            dset.time.gps.datetime,
        ],
        y_arrays=[dset.gdop, dset.pdop, dset.vdop, dset.hdop, dset.tdop],
        xlabel="Time [GPS]",
        ylabel="Dilution of precision",
        y_unit="",
        labels=["GDOP", "PDOP", "VDOP", "HDOP", "TDOP"],
        figure_path=figure_dir / f"plot_dop.{FIGURE_FORMAT}",
        opt_args={"figsize": (7, 4), "legend": True, "plot_to": "file"},
    )


def _plot_number_of_satellites(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Plot number of satellites

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """

    if "num_satellite_used" not in dset.fields:
        dset.add_float("num_satellite_used", val=gnss.get_number_of_satellites(dset))

    plot(
        x_arrays=[dset.time.gps.datetime, dset.time.gps.datetime],
        y_arrays=[dset.num_satellite_available, dset.num_satellite_used],
        xlabel="Time [GPS]",
        ylabel="Number of satellites",
        y_unit="",
        labels=["Available", "Used"],
        figure_path=figure_dir / f"plot_number_of_satellites.{FIGURE_FORMAT}",
        opt_args={"figsize": (7, 4), "legend": True, "marker": ",", "plot_to": "file", "plot_type": "plot"},
    )


def _plot_skyplot(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Plot skyplot

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """

    # Convert azimuth to range 0-360 degree
    azimuth = dset.site_pos.azimuth
    idx = azimuth < 0
    azimuth[idx] = 2 * np.pi + azimuth[idx]

    # Convert zenith distance from radian to degree
    zenith_distance = np.rad2deg(dset.site_pos.zenith_distance)

    # Generate x- and y-axis data per satellite
    x_arrays = []
    y_arrays = []
    labels = []
    for sat in dset.unique("satellite"):
        idx = dset.filter(satellite=sat)
        x_arrays.append(azimuth[idx])
        y_arrays.append(zenith_distance[idx])
        labels.append(sat)

    # Plot with polar projection
    # TODO: y-axis labels are overwritten after second array plot. Why? What to do?
    plot(
        x_arrays=x_arrays,
        y_arrays=y_arrays,
        xlabel="",
        ylabel="",
        y_unit="",
        labels=labels,
        figure_path=figure_dir / f"plot_skyplot.{FIGURE_FORMAT}",
        opt_args={
            "colormap": "tab20",
            "figsize": (7, 7.5),
            "legend": True,
            "legend_ncol": 6,
            "legend_location": "bottom",
            "plot_to": "file",
            "plot_type": "scatter",
            "projection": "polar",
            "title": "Skyplot\n Azimuth [deg] / Elevation[deg]",
            "xlim": [0, 2 * np.pi],
            "ylim": [0, 90],
            "yticks": (range(0, 90, 30)),  # sets 3 concentric circles
            "yticklabels": (map(str, range(90, 0, -30))),  # reverse labels from zenith distance to elevation
        },
    )


def _plot_satellite_elevation(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Plot satellite elevation

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """

    # Convert elevation from radian to degree
    elevation = np.rad2deg(dset.site_pos.elevation)

    # Limit x-axis range to rundate
    day_start, day_end = _get_day_limits(dset)

    # Generate x- and y-axis data per satellite
    x_arrays = []
    y_arrays = []
    labels = []

    for sat in dset.unique("satellite"):
        idx = dset.filter(satellite=sat)
        x_arrays.append(dset.time.gps.datetime[idx])
        y_arrays.append(elevation[idx])
        labels.append(sat)

    # Plot with scatter plot
    plot(
        x_arrays=x_arrays,
        y_arrays=y_arrays,
        xlabel="Time [GPS]",
        ylabel="Elevation [deg]",
        y_unit="",
        labels=labels,
        figure_path=figure_dir / f"plot_satellite_elevation.{FIGURE_FORMAT}",
        opt_args={
            "colormap": "tab20",
            "figsize": (7, 8),
            "legend": True,
            "legend_ncol": 6,
            "legend_location": "bottom",
            "plot_to": "file",
            "plot_type": "scatter",
            "title": "Satellite elevation",
            "xlim": [day_start, day_end],
        },
    )


def _plot_satellite_overview(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Plot satellite observation overview

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """
    figure_path = figure_dir / f"plot_satellite_overview.{FIGURE_FORMAT}"

    # Limit x-axis range to rundate
    day_start, day_end = _get_day_limits(dset)

    # Get time and satellite data from read and orbit stage
    time_read, satellite_read = _sort_by_satellite(
        _get_dataset(dset, stage="read", systems=dset.meta["obstypes"].keys())
    )
    time_orbit, satellite_orbit = _sort_by_satellite(
        _get_dataset(dset, stage="orbit", systems=dset.meta["obstypes"].keys())
    )
    time_edit, satellite_edit = _sort_by_satellite(
        _get_dataset(dset, stage="edit", systems=dset.meta["obstypes"].keys())
    )
    if not time_read:
        # NOTE: This is the case for concatencated Datasets, where "read" and "edit" stage data are not available.
        log.warn(f"No data for read stage available. Plot {figure_path} can not be plotted.")
        return 1

    # Generate plot
    plot(
        x_arrays=[time_read, time_orbit, time_edit],
        y_arrays=[satellite_read, satellite_orbit, satellite_edit],
        xlabel="Time [GPS]",
        ylabel="Satellite",
        y_unit="",
        # labels = ["Rejected in orbit stage", "Rejected in edit stage", "Kept observations"],
        colors=["red", "orange", "green"],
        figure_path=figure_path,
        opt_args={
            "colormap": "tab20",
            "figsize": (7, 6),
            "marker": "|",
            "plot_to": "file",
            "plot_type": "scatter",
            "title": "Overview over satellites",
            "xlim": [day_start, day_end],
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
    if dset_outlier == 1:
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
def _get_day_limits(dset: "Dataset") -> Tuple[datetime, datetime]:
    """Get start and end time for given run date

        Args:
            dset: A dataset containing the data.

        Returns:
            Start and end date. 
        """
    day_start = min(dset.time.datetime)
    day_end = max(dset.time.datetime)

    return day_start, day_end


def _get_outliers_dataset(dset: "Dataset") -> "Dataset":
    """Get dataset with outliers

    Args:
       dset:        A dataset containing the data.

    Returns:
       Dataset with outliers or status 1 if no data for "calculate" stage are available
    """

    # Get Dataset where no outliers are rejected
    dset_complete = dataset.Dataset(
        rundate=dset.analysis["rundate"],
        pipeline=dset.vars["pipeline"],
        stage="calculate",
        station=dset.vars["station"],
    )

    if dset_complete.num_obs == 0:
        # NOTE: This is the case for concatencated Datasets, where "calculate" stage data are not available.
        return 1

    # Get relative complement, which corresponds to "outlier" dataset
    dset_outliers = dset_complete.complement_with(dset, complement_by=["time", "satellite"])

    return dset_outliers


def _get_dataset(dset: "Dataset", stage: str, systems: Union[List[str], None] = None) -> "Dataset":
    """Get dataset for given stage

    Args:
       dset:        A dataset containing the data.
       systems:     List with GNSS identifiers (e.g. E, G, ...)

    Returns:
       Dataset for given stage
    """

    # Get Dataset
    # TODO: "label" should have a default value.
    dset_out = dset.read(rundate=dset.analysis["rundate"], pipeline=dset.vars["pipeline"], stage=stage, label="None")

    # Reject not defined GNSS observations
    if systems:
        systems = [systems] if isinstance(systems, str) else systems
        keep_idx = np.zeros(dset_out.num_obs, dtype=bool)
        for sys in systems:
            idx = dset_out.filter(system=sys)
            keep_idx[idx] = True
        dset_out.subset(keep_idx)

    return dset_out


def _sort_by_satellite(dset: "Dataset") -> Tuple[np.ndarray, np.ndarray]:
    """Sort time and satellite fields of dataset by satellite order

    Args: 
       dset:        A dataset containing the data.

    Returns:
        Tuple with ordered time array and satellite array
    """
    time = []
    satellite = []
    for sat in sorted(dset.unique("satellite"), reverse=True):
        idx = dset.filter(satellite=sat)
        time.extend(dset.time.gps.datetime[idx])
        satellite.extend(dset.satellite[idx])

    return time, satellite
