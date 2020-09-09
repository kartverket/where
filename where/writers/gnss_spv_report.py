"""Write report about a GNSS analysis run

Description:
------------

asdf


"""

# External library imports
import numpy as np

# Midgard imports
from where.data import position
from midgard.dev import plugins
from midgard.plot.matplotlib_extension import plot_scatter_subplots, plot

# Where imports
from where.lib import config
from where.writers._report import Report


FIGURE_DPI = 200
FIGURE_FORMAT = "png"


@plugins.register
def gnss_spv_report(dset: "Dataset") -> None:
    """Write report about a GNSS analysis run

    Args:
        dset:        A dataset containing the data.
    """

    # TODO: Better solution?
    if "station" not in dset.vars:  # necessary if called for example by ./where/tools/concatenate.py
        dset.vars["station"] = ""
        dset.vars["STATION"] = ""

    # Generate figure directory to save figures generated for GNSS report
    figure_dir = config.files.path("output_gnss_spv_report_figure", file_vars=dset.vars)
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    _plot_velocity(dset, figure_dir)
    # TODO_plot_residual(dset, figure_dir)
    _plot_dop(dset, figure_dir)

    # Generate GNSS report
    path = config.files.path(f"output_gnss_spv_report", file_vars=dset.vars)
    with config.files.open_path(path, create_dirs=True, mode="wt") as fid:
        rpt = Report(fid, rundate=dset.analysis["rundate"], path=path, description="GNSS SPV analysis")
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
    rpt.add_text("\n# GNSS site velocity analysis\n\n")

    # Site velocity
    rpt.add_figure(
        f"{figure_dir}/plot_timeseries_enu_vel.{FIGURE_FORMAT}",
        caption="Site velocity for East, North and Up direction",
        clearpage=True,
    )

    # Horizontal velocity
    rpt.add_figure(
        f"{figure_dir}/plot_horizontal_velocity.{FIGURE_FORMAT}", caption="Horizontal velocity", clearpage=True
    )

    # Plot 2D and 3D site velocity
    rpt.add_figure(
        f"{figure_dir}/plot_timeseries_2d_3d_vel.{FIGURE_FORMAT}",
        caption="Precision of dilution, 2D and 3D site velocity based on Doppler",
        clearpage=True,
    )

    #    #
    #    # Residual
    #    #
    #    rpt.add_text("\n# GNSS residual\n\n")
    #
    #    # Plot residuals
    #    rpt.add_figure(
    #            f"{figure_dir}/plot_residual.{FIGURE_FORMAT}",
    #            caption="Post-fit residuals.",
    #            clearpage=True
    #    )

    #
    # Dilution of precision (DOP)
    #
    rpt.add_text("\n# Dilution of precision\n\n")

    # Plot DOP
    rpt.add_figure(f"{figure_dir}/plot_dop.{FIGURE_FORMAT}", caption="Dilution of precision.", clearpage=True)


#
# PLOT FUNCTIONS
#
def _plot_velocity(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Plot site position plots

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """

    plot_scatter_subplots(
        x_array=dset.time.gps.datetime,
        y_arrays=[dset.e_vel, dset.n_vel, dset.u_vel],
        xlabel="Time [GPS]",
        ylabels=["East", "North", "Up"],
        colors=["steelblue", "darkorange", "limegreen"],
        y_units=["m/s", "m/s", "m/s"],
        figure_path=figure_dir / f"plot_timeseries_enu_vel.{FIGURE_FORMAT}",
        opt_args={
            "figsize": (6, 6),
            "plot_to": "file",
            "sharey": True,
            "title": "Site velocity",
            "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
        },
    )

    plot_scatter_subplots(
        x_array=dset.e_vel,
        y_arrays=[dset.n_vel],
        xlabel="East [m/s]",
        ylabels=["North"],
        y_units=["m/s"],
        figure_path=figure_dir / f"plot_horizontal_velocity.{FIGURE_FORMAT}",
        opt_args={
            "grid": True,
            "figsize": (6, 6),
            "histogram": "x, y",
            "plot_to": "file",
            "title": "Horizontal velocity",
            "xlim": [-0.2, 0.2],
            "ylim": [-0.2, 0.2],
        },
    )

    plot_scatter_subplots(
        x_array=dset.time.gps.datetime,
        y_arrays=[dset.pdop, dset["vel_2d"], dset["vel_3d"]],
        xlabel="Time [GPS]",
        ylabels=["PDOP", "2D", "3D"],
        colors=["steelblue", "darkorange", "limegreen"],
        y_units=[None, "m/s", "m/s"],
        figure_path=figure_dir / f"plot_timeseries_2d_3d_vel.{FIGURE_FORMAT}",
        opt_args={
            "figsize": (6, 6),
            "plot_to": "file",
            "sharey": False,
            "title": "2D and 3D site velocity",
            "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
        },
    )


def _plot_residual(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Plot residual plot

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """
    plot(
        x_arrays=dset.time.gps.datetime,
        y_arrays=dset.residual,
        xlabel="Time [GPS]",
        ylabel="Post-fit residual",
        y_unit="meter",
        colors=["dodgerblue"],
        figure_path=figure_dir / f"plot_residual.{FIGURE_FORMAT}",
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
        x_arrays=[dset.time.gps.datetime, dset.time.gps.datetime, dset.time.gps.datetime, dset.time.gps.datetime],
        y_arrays=[dset.gdop, dset.pdop, dset.hdop, dset.vdop],
        xlabel="Time [GPS]",
        ylabel="Dilution of precision",
        y_unit="meter",
        labels=["GDOP", "PDOP", "HDOP", "VDOP"],
        figure_path=figure_dir / f"plot_dop.{FIGURE_FORMAT}",
        opt_args={"figsize": (7, 4), "legend": True, "plot_to": "file"},
    )
