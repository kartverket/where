"""Generate SISRE plots

Description:
------------

Plots can be configured via configuration file where_pipeline_sisre.conf. See section 'sisre_plot'.

"""
# Standard library imports
from collections import namedtuple
from pathlib import PosixPath
import textwrap
from typing import Any, Dict, List, Union

# External library imports
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


# Midgard imports
from midgard.dev import plugins
from midgard.plot.matplotext import MatPlotExt
from midgard.math.unit import Unit

# Where imports
from where.lib import config
from where.lib import log
from where.lib import util


PlotConfig = namedtuple("PlotConfig", ["label", "statistic", "format"])
PlotConfig.__doc__ = """A convenience class for defining a field for subplot configuration

    Args:
        label (str):      Label
        statistic(bool):  Flag for plotting statistical information in title or not
        format (str):     Format of statistical values (e.g. mean or rms) in title
    """

SubplotConfig = namedtuple("SubplotConfig", ["ylabel", "color", "ydata"])
SubplotConfig.__doc__ = """A convenience class for defining a field for subplot configuration

    Args:
        ylabel (str):           Y-axis label
        color (str):            Color of scatter plot
        ydata (numpy.ndarray):  Y-axis data
    """


# Define configuration of plots
FIGURE_FORMAT = "png"
GNSS_NAME = {"C": "BeiDou", "E": "Galileo", "G": "GPS", "I": "IRNSS", "J": "QZSS", "R": "GLONASS"}

# TODO: How to get statistical format specification to work?
PLOTCONFIG = {
    "age_of_ephemeris": PlotConfig("Age of ephemeris", True, ".0f"),
    "sisre": PlotConfig("SISE", True, ".2f"),
    "sisre_orb": PlotConfig("orbit-only SISE", True, ".2f"),
    "orb_diff_3d": PlotConfig("3D orbit error", True, ".2f"),
    "clk_diff": PlotConfig("Clock difference \n without mean", True, ".2f"),
    "clk_diff_with_dt_mean": PlotConfig("Clock difference", True, ".2f"),
    "bias_brdc": PlotConfig("Satellite bias of broadcast clocks", True, ".0f"),
    "bias_precise": PlotConfig("Satellite bias of precise clocks", True, ".0f"),
    "used_iode": PlotConfig("IODE", False, ".0f"),
}

UNIT_TITLE = {"": "", "meter": "m", "minute": "min", "second": "s"}
UNIT_YLABEL = {"": "", "meter": "m", "minute": "min", "second": "s"}


@plugins.register
def sisre_plot(dset):
    """Write SISRE report

    Args:
        dset (Dataset):       A dataset containing the data.
    """

    # Satellites and systems to plot
    satellites = config.where.sisre_plot.get("satellites", default=dset.unique("satellite")).list
    satellites = satellites if satellites else dset.unique("satellite")
    systems = config.where.sisre_plot.get("systems", default="").list
    systems = systems if systems else dset.unique("system")

    # Overwrite configuration with command options
    satellite = util.read_option_value("--satellite", default=None)
    satellites = [satellite] if satellite else satellites
    title = util.read_option_value("--title", default=None)

    # Set matplotlib configuration
    options = _set_plot_config(title)

    # Get fields to plot
    fields = config.where.sisre_plot.get("fields", default="").list
    if not fields:
        log.fatal(f"No fields to plot.")

    if options["subplot"]:
        _plot_scatter_subplots(dset, fields, systems, satellites, options)

    for field in fields:
        if field == "age_of_ephemeris":
            dset[field][:] = dset[field] * Unit.second2minute
            dset.set_unit(field, "minute")

        _plot_scatter_field(dset, field, systems, satellites, options)


def _set_plot_config(title: str) -> Dict[str, Any]:
    """Set matplotlib configuration by reading sisre_plot configuration

    Args:
        title: Title given via option. Overwrites configuration file definition.
    """
    # Change fontsize of labels
    fontsize = config.where.sisre_plot.get("fontsize", default=10).str
    fontsize = 9 if not fontsize else int(fontsize)
    plt.rcParams["axes.titlesize"] = fontsize
    plt.rcParams["axes.labelsize"] = fontsize
    plt.rcParams["font.size"] = fontsize
    plt.rcParams["figure.titlesize"] = fontsize
    plt.rcParams["legend.fontsize"] = fontsize
    plt.rcParams["xtick.labelsize"] = fontsize
    plt.rcParams["ytick.labelsize"] = fontsize

    # Define additional matplotlib/plot configuration
    options = {
        "alpha": config.where.sisre_plot.get("alpha", default=1).int,
        "dpi": config.where.sisre_plot.get("dpi", default=200).int,
        "color": config.where.sisre_plot.get("color", default="").str,
        "colormap": config.where.sisre_plot.get("colormap", default="tab10").str,
        "figsize": tuple([int(v) for v in config.where.sisre_plot.get("figsize", default="6,4").list]),
        "figsize_subplot": tuple([int(v) for v in config.where.sisre_plot.get("figsize_subplot", default="7,5").list]),
        "fontsize": fontsize,
        "fsize_subtitle": config.where.sisre_plot.get("fontsize_subtitle", default=7).int,
        "legend": config.where.sisre_plot.get("legend", default=True).bool,
        "marker": config.where.sisre_plot.get("marker", default=".").str,
        "markersize": config.where.sisre_plot.get("markersize", default=9).int,
        "subplot": config.where.sisre_plot.get("subplot", default=True).bool,
        "title": title if title else config.where.sisre_plot.get("title", default=" ").str,
    }

    return options


def _get_figure_path(
    dset: "Dataset", field: str, sys: str, satellite: Union[List[str], str] = "", subplot: bool = False
) -> PosixPath:
    """Get figure path and generate figure directory
    Args:
       dset:       A dataset containing the data.
       field:      Dataset field.
       sys:        GNSS identifier.
       subplot:    Plot subplot or not. 

    Returns:
       Figure path
    """
    if len(satellite) == 1:
        satellite = satellite[0]
    else:
        satellite = ""
    dset.vars.update(
        {"field": field, "format": FIGURE_FORMAT, "system": GNSS_NAME[sys].lower(), "satellite": satellite.lower()}
    )
    file_key = "output_sisre_subplot" if subplot else "output_sisre_plot"
    figure_path = config.files.path(file_key, file_vars={**dset.vars, **dset.analysis})
    if not figure_path.parent.exists():
        figure_path.parent.mkdir(parents=True, exist_ok=True)
    log.info(f"Plot figure: {figure_path}")
    return figure_path


def _get_statistic_text(data: np.ndarray, unit: str = "", width: int = 150) -> str:
    """Get statistical text string

    Args:
        data:   Array with data.
        unit:   Unit of data
        width:  Maximal width of statistical text 


    Return:
        text representing statistical information      
    """
    mpe = MatPlotExt()
    stat_text = ", ".join(mpe.get_statistic(data, unit=unit))
    return textwrap.fill(stat_text, width=width, break_long_words=False, break_on_hyphens=True)


def _insert_kartverket_logo(fig):
    # TODO: Does not work!!
    logo_path = "/home/dahmic/texmf/tex/figure/kartverket_staende.png"
    img = plt.imread(logo_path)
    fig.figimage(img, xo=100, yo=100, alpha=0.15, origin="upper", zorder=1)  # , extent=[0,0,1,1])


#
# SCATTER SUPPLOTS
#
def _plot_scatter_subplots(
    dset: "Dataset", fields: List[str], systems: List[str], satellites: List[str], options: Dict[str, Any]
) -> None:
    """Generate scatter subplot
    Args:
       dset:       A dataset containing the data.
       field:      Dataset field.
       options:    Dictionary with options, which overwrite default plot configuration.
    """
    # Generate plot for each GNSS
    for sys in systems:
        idx = dset.filter(system=sys)
        xlim_min = None
        xlim_max = None

        # Define colormap
        colormap = getattr(cm, options["colormap"])
        colors = colormap(np.linspace(0, 1, len(dset.unique("satellite", idx=idx))))

        fig, axes = plt.subplots(len(fields), 1, sharex=True, sharey=False, figsize=options["figsize_subplot"])
        fig.suptitle(f"{options['title']}", y=1.0)

        for ax, field in zip(axes, fields):

            if field == "age_of_ephemeris":
                dset[field][:] = dset[field] * Unit.second2minute
                dset.set_unit(field, "minute")

            unit = "" if dset.unit(field) is None else dset.unit(field)[0]

            # Plot all satellites or only defined ones
            for sat, color in zip(sorted(dset.unique("satellite", idx=idx)), colors):
                if sat not in satellites:
                    continue
                idx_sat = dset.filter(satellite=sat)

                # Note: Normally for each satellite the correct x-limit range is chosen, but it is not the case that
                #      the range is the same for all satellites.
                if xlim_min is None:
                    xlim_min = min(dset.time.gps.datetime[idx_sat])
                    xlim_max = max(dset.time.gps.datetime[idx_sat])
                else:
                    xlim_min = (
                        min(dset.time.gps.datetime[idx_sat])
                        if xlim_min > min(dset.time.gps.datetime[idx_sat])
                        else xlim_min
                    )
                    xlim_max = (
                        max(dset.time.gps.datetime[idx_sat])
                        if xlim_max < max(dset.time.gps.datetime[idx_sat])
                        else xlim_max
                    )
                options.update({"xlim": [xlim_min, xlim_max]})

                # Overwrite colormap color selection
                color = options["color"] if options["color"] else color

                if np.sum(dset[field][idx_sat]) != 0:
                    mpe = MatPlotExt()
                    mpe.plot_subplot_row(
                        ax,
                        x_array=dset.time.gps.datetime[idx_sat],
                        y_array=dset[field][idx_sat],
                        ylabel=f"{PLOTCONFIG[field].label}",
                        y_unit=UNIT_YLABEL[unit],
                        label=sat,
                        color=color,
                        options=options,
                    )

            # Overwrite subtitle with statistical information over all satellites
            statistic = _get_statistic_text(dset[field], unit=unit)
            ax.set_title(statistic, fontsize=options["fsize_subtitle"], horizontalalignment="center")

        # Plot x-axis label only once below the last subplot row
        ax.set(xlabel="Time [GPS]")

        # Rotates and right aligns the x labels, and moves the bottom of the axes up to make room for them
        fig.autofmt_xdate()

        # Plot legend on the right side of the figure by getting labels only from the last subplot row
        if options["legend"]:
            handles, labels = ax.get_legend_handles_labels()
            fig.legend(
                handles, 
                labels, 
                loc="upper center", 
                ncol=7, 
                bbox_to_anchor=(0.54, -0.005),
                fontsize=8,
                fancybox=True, 
                shadow=True,
            )

        # Adjust subplot axes (to place legend on the right side of the figure)
        #fig.subplots_adjust(right=0.89, top=0.91)

        # Automatically adjusts subplot parameteres so that the subplot(s) fits in to the figure area
        # Note: Legend is not taken into account.
        fig.tight_layout()

        # Save figure
        plt.savefig(
            _get_figure_path(dset, field, sys, satellite=satellites, subplot=True),
            dpi=options["dpi"],
            bbox_inches="tight",
        )

        # Clear the current figure
        plt.clf()


def _plot_scatter_field(
    dset: "Dataset", field: List[str], systems: List[str], satellites: List[str], options: Dict[str, Any]
) -> None:
    """Scatter plot of given field

    Args:
       dset (Dataset):           A dataset containing the data.
       field (str):              Dataset field.
       legend(bool):             Plot legend or not.

    Returns:
        None
    """
    unit = "" if dset.unit(field) is None else dset.unit(field)[0]

    # Generate plot for each GNSS
    for sys in systems:
        idx = dset.filter(system=sys)

        # Define colormap
        colormap = getattr(cm, options["colormap"])
        colors = colormap(np.linspace(0, 1, len(dset.unique("satellite", idx=idx))))

        fig = plt.figure(figsize=options["figsize"])
        # TODO: does not work _insert_kartverket_logo(fig)

        # Plot all satellites or only defined ones
        for sat, color in zip(sorted(dset.unique("satellite", idx=idx)), colors):
            if sat not in satellites:
                continue
            idx_sat = dset.filter(system=sys, satellite=sat)

            # Overwrite colormap color selection
            color = options["color"] if options["color"] else color

            if np.sum(dset[field][idx_sat]) != 0:
                plt.scatter(
                    list(dset.time.gps.datetime[idx_sat]),
                    list(dset[field][idx_sat]),
                    color=color,
                    label=sat,
                    alpha=options["alpha"],
                    marker=options["marker"],
                    s=options["markersize"],
                )
            sat_statistic = _get_statistic_text(dset[field][idx_sat], unit=unit)

        if options["legend"]:
            # plt.legend(bbox_to_anchor=(1.2, 1), ncol=1)
            #plt.legend(
            #    bbox_to_anchor=(1.04, 1), 
            #    borderaxespad=0.0, 
            #    loc="upper left", 
            #    ncol=1,
            #    fancybox=True, 
            #    shadow=True,
            #)

            plt.legend(
                loc="upper center", 
                ncol=6, 
                bbox_to_anchor=(0.50, -0.20), 
                fontsize=8,
                fancybox=True, 
                shadow=True,
            )

        ylabel = f"{PLOTCONFIG[field].label} [{UNIT_YLABEL[unit]}]" if UNIT_YLABEL[unit] else f"{PLOTCONFIG[field].label}"
        plt.ylabel(ylabel)
        plt.xlim([min(dset.time.gps.datetime[idx]), max(dset.time.gps.datetime[idx])])
        plt.xlabel("Time [GPS]")

        # Plot statistic only if one satellite or all satellites are given
        if len(satellites) == 1:
            statistic = sat_statistic
        elif len(satellites) == len(dset.unique("satellite")):
            statistic = _get_statistic_text(dset[field][idx], unit=unit)
        else:
            statistic = ""

        # Plot title
        if PLOTCONFIG[field].statistic:
            plt.title(options["title"], y=1.05)
            plt.suptitle(statistic, y=0.93, fontsize=options["fsize_subtitle"])
        else:
            plt.title(options["title"])

        # Rotates and right aligns the x labels, and moves the bottom of the axes up to make room for them
        fig.autofmt_xdate()

        # Automatically adjusts subplot parameteres so that the subplot(s) fits in to the figure area
        # Note: Legend is not taken into account.
        fig.tight_layout()

        # Adjust plot axes (to place title and legend correctly)
        fig.subplots_adjust(right=0.89, top=0.90)

        # Save figure
        plt.savefig(_get_figure_path(dset, field, sys, satellite=satellites), dpi=options["dpi"], bbox_inches="tight")

        # Clear the current figure
        plt.clf()
