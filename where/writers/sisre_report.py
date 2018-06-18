"""Write report about a SISRE analysis run (only for one session and not several sessions (stations))

Description:
------------

asdf




"""
# Standard library imports
from collections import namedtuple
import getpass
from datetime import datetime
import re
import textwrap

# External library imports
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# WHERE imports
import where
from where import apriori
from where import cleaners
from where.lib import config
from where.lib import enums
from where.lib import files
from where.lib import gnss
from where.lib import plugins

gnss_name = {"C": "BeiDou", "E": "Galileo", "G": "GPS", "I": "IRNSS", "J": "QZSS", "R": "GLONASS"}

SubplotConfig = namedtuple("SubplotConfig", ["ylabel", "color", "ydata"])
SubplotConfig.__doc__ = """A convenience class for defining a field for subplot configuration

    Args:
        ylabel (str):           Y-axis label
        color (str):            Color of scatter plot
        ydata (numpy.ndarray):  Y-axis data
    """


AxhlineConfig = namedtuple("AxhlineConfig", ["type", "y_value", "color"])
AxhlineConfig.__doc__ = """A convenience class for defining matplotlib axhline configuration

    Args:
        type (str):             Type (e.g. GNSS identifier E or G for Galileo or GPS)
        y_value (float):        Y-value for horizontal line to be plotted in [m]
        color (str):            Color of horizontal line
    """


@plugins.register
def sisre_report(dset):
    """Write SISRE report

    Args:
        dset:       Dataset, a dataset containing the data.
    """
    write_level = config.tech.get("write_level", default="operational").as_enum("write_level")

    with files.open(file_key="output_gnss_session_report", file_vars=dset.vars, mode="wt") as fid:
        _header(fid)
        _write_config(fid)
        fid.write("\n# Satellite status\n\n")
        _unhealthy_satellites(fid, dset)
        _eclipse_satellites(fid, dset)

        # Generate figure directory to save figures generated for SISRE report
        figure_dir = files.path("output_sisre_report_figure", file_vars=dset.vars)
        figure_dir.mkdir(parents=True, exist_ok=True)
        if write_level <= enums.get_value("write_level", "detail"):
            _plot_scatter_satellite_bias(fid, figure_dir, dset)
            _plot_scatter_orbit_and_clock_differences(fid, figure_dir, dset)
        _plot_histogram_sisre(fid, figure_dir, dset)
        _plot_scatter_sisre(fid, figure_dir, dset)
        _statistics(fid, figure_dir, dset)


def _eclipse_satellites(fid, dset):
    """Write overview over satellites in eclipse

    Args:
       fid:
       dset:
    """
    # TODO: time period of eclipting satellites needed
    brdc = apriori.get(
        "orbit",
        rundate=dset.rundate,
        time=dset.time,
        satellite=tuple(dset.satellite),
        system=tuple(dset.system),
        station=dset.dataset_name.upper(),
        apriori_orbit="broadcast",
    )

    fid.write(
        "{:25s} = {sats:47s}\n"
        "".format("Satellites in eclipse", sats=" ".join(gnss.check_satellite_eclipse(brdc.dset)))
    )


def _header(fid):
    title = "Where"

    fid.write(
        """---
title: {title}
author: Where v{version} [{user}]
date: {nowdate:%Y-%m-%d}
---
""".format(
            # title=title["text"], version=where.__version__, user=config.analysis.user.str, nowdate=datetime.now()
            title=title,
            version=where.__version__,
            user=getpass.getuser(),
            nowdate=datetime.now(),  # TODO: Better solution?
        )
    )

    fid.write("\\newpage\n")
    fid.write(
        """#SISRE analysis\n

For the SISRE analysis it is common to apply the average contribution over all points of the Earth within the visibility cone of the satellite (Montenbruck el al., 2014), which is called global averaged SISRE. The SISRE analysis in Where is based on the global averaged SISRE:

\\begin{equation}   
     \\text{SISRE}^s = \sqrt{(w_r \cdot \Delta r^s - \Delta t^s)^2 + w_{a,c}^2 \cdot (\Delta {a^s}^2 + \Delta {c^s}^2)}
  \\label{eq:sisre}  
\\end{equation}

\\begin{equation}   
     \\text{SISRE(orb)}^s = \sqrt{w_r^2 \cdot \Delta {r^s}^2 + w_{a,c}^2 \cdot (\Delta {a^s}^2 + \Delta {c^s}^2)}
  \\label{eq:sisre_orb}  
\\end{equation}

\\begin{tabular}{lll}
with\\\\
 & $\\text{SISRE}^s$        & Global averaged SISRE for satellite $s$, \\\\
 & $\\text{SISRE(orb)}^s$   & Orbit-only SISRE for satellite $s$, \\\\
 &$\Delta a^s$, $\Delta c^s$, $\Delta r^s$   & Satellite coordinate differences between broadcast and precise\\\\
 &                                           & ephemeris in along-track, cross-track and radial for satellite $s$,\\\\
 &$\Delta t^s$              & Satellite clock correction difference related to CoM of \\\\
 &                          & satellite $s$ and corrected for satellite biases and averaged \\\\
 &                          & clock offset in each epoch,\\\\
 &$w_r$                     & SISRE weight factor for radial errors (see Table 1),\\\\
 &$w_{a,c}$                 & SISRE weight factor for along-track and cross-track errors\\\\
 &                          & (see Table 1).\\\\
\\end{tabular}

It should be noted that we have neglected the uncertainty of precise ephemeris and clocks by the determination of SISRE.

The SISRE analysis is carried out on daily basis. Each daily solution is cleaned for outliers. The outliers are detected and rejected for each day iteratively using a 3-sigma threshold determined for each satellite. After each iteration the SISRE results are again recomputed.

\\begin{table}[!ht]
\\begin{center}
  \\begin{tabular}[c]{lll}
    \\hline
    System       & $w_r$  & $w_{a,c}$ \\\\
    \\hline
    GPS          & $0.980$ & $0.139$ \\\\
    Galileo      & $0.984$ & $0.124$ \\\\
    \\hline
  \\end{tabular}
  \\caption[The global averaged SISRE weight factors for radial ($w_r$) and along-track and cross-track errors ($w_{a,c}$). The weight factors are given for an elevation mask of $5^{\circ}$ based on Table 4 in \cite{montenbruck2018}.]{The global averaged SISRE weight factors for radial ($w_r$) and along-track and cross-track errors ($w_{a,c}$). The weight factors are given for an elevation mask of $5^{\circ}$ based on Table 4 in Montenbruck et al. 2018.}
  \\label{tab:sisre_weight_factors}
\\end{center}
\\end{table}

\\newpage\n """
    )


def _get_satellite_type(dset, satellite):
    idx = dset.filter(satellite=satellite)
    return set(dset.satellite_type[idx]).pop()


def _plot_bar_dataframe_columns(fid, figure_dir, df, field, extra_row_names, column="rms"):
    """Generate bar plot of given dataframe columns (colored and ordered by satellite type)
    """
    fontsize = 12
    df_reduced = df.drop(extra_row_names)  # Remove extra rows

    # Assign to each satellite type a color
    colors = dict()
    # TODO: Better handling of color definition?
    # color_def = ['cornflowerblue', 'firebrick', 'violet', 'gold', 'limegreen', 'deepskyblue', 'orangered']
    color_def = ["red", "tomato", "lightsalmon", "blue", "royalblue", "deepskyblue", "paleturquoise"]
    # color_def = ['C'+str(idx) for idx in range(0,10)]
    for type_ in sorted(set(df_reduced.type)):
        colors.update({type_: color_def.pop()})

    # Generate bar plot
    df_color = df_reduced["type"].apply(lambda x: colors[x])
    fig_width = len(df_reduced.index) / 4 if len(df_reduced.index) > 30 else 6.4
    ax = df_reduced[column].plot(kind="bar", color=df_color, width=0.8, figsize=(fig_width, fig_width / 1.33))
    ax.set_xlabel("Satellite", fontsize=fontsize)
    ax.set_ylabel(f"{field.upper()} {column.upper()} [m]", fontsize=fontsize)

    # Make legend
    satellite_type_patch = [mpatches.Patch(color=v, label=k) for k, v in sorted(colors.items())]
    ax.legend(handles=satellite_type_patch, bbox_to_anchor=(1.04, 1), loc=2, borderaxespad=0., ncol=1)

    plt.tight_layout()
    plt.savefig(figure_dir / f"plot_bar_{field}_{column}.pdf")
    plt.clf()  # clear the current figure

    fid.write(
        f"![{field.upper()} {column.upper()} for all satellites sorted by satellite type]({figure_dir}/plot_bar_{field}_{column}.pdf)\n"
    )
    fid.write("\\newpage\n")


#
# BAR STACKED PLOTS
#
def _plot_bar_stacked(fid, df_sisre, df_sisre_orb, figure_path, xlabel="", xticks_rotation=None, axhline=None):
    """Generate bar plot of given dataframe columns (colored and ordered by satellite type)
    """
    fontsize = 12

    # Generate new dataframe with columns 'orbit-only SISRE' and 'clock-only SISRE'
    data = np.hstack(
        (
            np.array([df_sisre_orb.rms]).T,
            (np.array([df_sisre.rms]) - np.array([df_sisre_orb.rms])).T,
            (np.array([df_sisre.percentile]) - np.array([df_sisre.rms])).T,
        )
    )
    df_merged = pd.DataFrame(
        data=data,
        index=df_sisre.index,
        columns=["orbit-only SISRE RMS", "clock-only SISRE RMS", "95th percentile SISRE"],
    )

    # Generate bar plot
    fig_width = len(df_merged.index) / 4 if len(df_merged.index) > 30 else 6.4
    ax = df_merged.plot(kind="bar", stacked=True, width=0.8, figsize=(fig_width, fig_width / 1.33))
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(f"Accuracy [m]", fontsize=fontsize)
    if axhline is not None:
        for idx in range(0, len(axhline)):
            plt.axhline(y=axhline[idx].y_value, linewidth=2, color=axhline[idx].color)
    # ax.legend(bbox_to_anchor=(1.04, 1), loc=2, borderaxespad=0., ncol=1)
    ax.legend(bbox_to_anchor=(1, 1.15), loc=1, borderaxespad=0., ncol=2)

    if xticks_rotation is not None:
        plt.xticks(rotation=xticks_rotation)
    plt.tight_layout()
    plt.savefig(figure_path)
    plt.clf()  # clear the current figure


def _plot_bar_stacked_sisre(fid, field_dfs, extra_row_names, figure_dir):

    figure_path_threshold = figure_dir / f"plot_bar_stacked_sisre_extra_rows_threshold.pdf"
    figure_path = figure_dir / f"plot_bar_stacked_sisre_extra_rows.pdf"

    # Define SISRE thresholds
    sisre_threshold = (AxhlineConfig("E", 2, "red"),)  # TODO: Handling of several GNSS

    # Remove satellite rows
    df_sisre_extra = field_dfs["sisre"].drop(set(field_dfs["sisre"].index) - set(extra_row_names))
    df_sisre_orb_extra = field_dfs["sisre_orb"].drop(set(field_dfs["sisre_orb"].index) - set(extra_row_names))

    _plot_bar_stacked(
        fid,
        df_sisre_extra,
        df_sisre_orb_extra,
        figure_path=figure_path_threshold,
        xticks_rotation=20,
        axhline=sisre_threshold,
    )
    _plot_bar_stacked(
        fid,
        df_sisre_extra,
        df_sisre_orb_extra,
        figure_path=figure_path_threshold,
        xticks_rotation=20,
        axhline=sisre_threshold,
    )
    _plot_bar_stacked(fid, df_sisre_extra, df_sisre_orb_extra, figure_path=figure_path, xticks_rotation=20)

    fid.write(
        f"![Blue bars indicate orbit-only SISRE RMS, orange bars clock-only SISRE RMS and green bars 95th percentile SISRE]({figure_path_threshold})\n"
    )
    fid.write(
        f"![Blue bars indicate orbit-only SISRE RMS, orange bars clock-only SISRE RMS and green bars 95th percentile SISRE]({figure_path})\n"
    )
    fid.write("\\newpage\n")


def _plot_bar_stacked_sisre_satellites(fid, field_dfs, extra_row_names, figure_dir):

    figure_path = figure_dir / f"plot_bar_stacked_sisre.pdf"

    # Remove extra rows
    df_sisre = field_dfs["sisre"].drop(extra_row_names)
    df_sisre_orb = field_dfs["sisre_orb"].drop(extra_row_names)

    _plot_bar_stacked(fid, df_sisre, df_sisre_orb, figure_path=figure_path, xlabel="Satellite")

    fid.write(
        f"![Blue bars indicate orbit-only SISRE RMS, orange bars clock-only SISRE RMS and green bars 95th percentile SISRE]({figure_path})\n"
    )
    fid.write("\\newpage\n")


#
# HISTOGRAM PLOTS
#
def _plot_histogram_subplot(data, axis, system):
    axis.hist(data, normed=True, bins=30)
    axis.set(xlabel="SISRE [m]", ylabel="Frequency")
    axis.set_title(f"{gnss_name[system]}")
    mean = np.mean(data)
    std = np.std(data)
    axis.text(
        0.98,
        0.98,
        f"$mean={mean:5.3f}\ \pm {std:5.3f}$ m",
        horizontalalignment="right",
        verticalalignment="top",
        transform=axis.transAxes,
    )


def _plot_histogram_sisre(fid, figure_dir, dset):
    import math

    # TODO: How to handle for example 3 subplots?
    nrows = math.ceil(len(dset.unique("system")) / 2)
    ncols = 2 if len(dset.unique("system")) > 1 else 1
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, squeeze=False)
    for sys, ax in zip(dset.unique("system"), axes.flatten()):
        idx = dset.filter(system=sys)
        _plot_histogram_subplot(dset.sisre[idx], ax, sys)
    plt.savefig(figure_dir / f"plot_histogram_sisre.pdf")
    plt.clf()  # clear the current figure

    fid.write(f"![Histrogram of SISRE results]({figure_dir}/plot_histogram_sisre.pdf)\n")
    fid.write("\\newpage\n")


#
# SCATTER PLOTS
#
def _plot_scatter_satellite_bias(fid, figure_dir, dset):

    for sys in dset.unique("system"):
        idx = dset.filter(system=sys)

        for field, orbit in {"bias_brdc": "Broadcast", "bias_precise": "Precise"}.items():
            if np.sum(dset[field][idx]) != 0:
                plt.scatter(dset.time.gps.datetime[idx], dset[field][idx])
                plt.ylabel(f"{orbit} satellite bias [m]")
                plt.xlim([min(dset.time.gps.datetime[idx]), max(dset.time.gps.datetime[idx])])
                plt.xlabel("Time [GPS]")
                plt.title(f"{gnss_name[sys]}")
                plt.savefig(figure_dir / f"plot_scatter_{field}_{gnss_name[sys].lower()}.pdf")
                plt.clf()  # clear the current figure

                fid.write(
                    f"![{orbit} satellite bias]({figure_dir}/plot_scatter_{field}_{gnss_name[sys].lower()}.pdf)\n"
                )
                fid.write("\\newpage\n")


#
# SCATTER SUPPLOTS
#
def _plot_scatter_subplots(xdata, subplots, figure_path, xlabel="", title=""):
    marker = "."  # point marker type

    fig, axes = plt.subplots(len(subplots), 1, sharex=True, sharey=True)
    fig.set_figheight(9)  # inches
    fig.suptitle(f"{title}", y=1.0)
    for idx, ax in enumerate(axes):
        ax.set(ylabel=subplots[idx].ylabel)
        ax.set_xlim([min(xdata), max(xdata)])  # otherwise time scale of x-axis is not correct -> Why?
        text = f"mean $ = {np.mean(subplots[idx].ydata):.2f} \pm {np.std(subplots[idx].ydata):.2f}$ m"
        ax.text(0.98, 0.98, text, horizontalalignment="right", verticalalignment="top", transform=ax.transAxes)
        ax.scatter(xdata, subplots[idx].ydata, marker=marker, color=subplots[idx].color)
        ax.set(xlabel=xlabel)

    fig.autofmt_xdate()  # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    plt.tight_layout()
    plt.savefig(figure_path)
    plt.clf()  # clear the current figure


def _plot_scatter_orbit_and_clock_differences(fid, figure_dir, dset):

    for sys in dset.unique("system"):
        idx = dset.filter(system=sys)

        # Define configuration of subplots
        subplots = (
            SubplotConfig("Δalong-track [m]", "paleturquoise", dset.orb_diff.itrs[:, 0][idx]),
            SubplotConfig("Δcross-track [m]", "deepskyblue", dset.orb_diff.itrs[:, 1][idx]),
            SubplotConfig("Δradial [m]", "royalblue", dset.orb_diff.itrs[:, 2][idx]),
            SubplotConfig("Δclock [m]", "tomato", dset.clk_diff_sys[idx]),
        )

        _plot_scatter_subplots(
            dset.time.gps.datetime[idx],
            subplots,
            figure_path=figure_dir / f"plot_scatter_orbit_clock_differences_{gnss_name[sys].lower()}.pdf",
            xlabel="Time [GPS]",
            title=f"{gnss_name[sys]}",
        )

        fid.write(
            f"![Broadcast ephemeris errors for along-track, cross-track and radial orbit errors and clock errors]({figure_dir}/plot_scatter_orbit_clock_differences_{gnss_name[sys].lower()}.pdf)\n"
        )
        fid.write("\\newpage\n")


def _plot_scatter_sisre(fid, figure_dir, dset):

    for sys in dset.unique("system"):
        idx = dset.filter(system=sys)

        # Define configuration of subplots
        subplots = (
            SubplotConfig("orbit-only SISRE [m]", "paleturquoise", dset.sisre_orb[idx]),
            SubplotConfig("clock-only SISRE [m]", "deepskyblue", dset.sisre[idx] - dset.sisre_orb[idx]),
            SubplotConfig("SISRE [m]", "royalblue", dset.sisre[idx]),
        )

        _plot_scatter_subplots(
            dset.time.gps.datetime[idx],
            subplots,
            figure_path=figure_dir / f"plot_scatter_sisre_{gnss_name[sys].lower()}.pdf",
            xlabel="Time [GPS]",
            title=f"{gnss_name[sys]}",
        )

        fid.write(
            f"![Orbit-only SISRE, clock-only SISRE and SISRE]({figure_dir}/plot_scatter_sisre_{gnss_name[sys].lower()}.pdf)\n"
        )
        fid.write("\\newpage\n")


def _statistics_satellite(fid, figure_dir, dset):
    """Generate statistics for each field and satellite

    Args:
        fid (_io.TextIOWrapper):         File object
        figure_dir (pathlib.PosixPath):  Figure directory path
        dset (Dataset):                  Dataset
    """
    fields = ["sisre", "sisre_orb", "orb_diff_3d", "clk_diff_sys"]
    rms = lambda x: np.sqrt(np.mean(np.square(x)))
    percentile = lambda x: np.percentile(x, 95)
    functions = [rms, np.mean, np.std, np.min, np.max, percentile]
    columns = ["type", "rms", "mean", "std", "min", "max", "percentile"]
    field_dfs = dict()

    # Generate field DataFrames with the satellites as indices and functional values (rms, mean, ...) as columns
    #
    # Example:
    #           type       rms      mean       std       min       max  percentile
    # E11  GALILEO-1  0.175557  0.138370  0.108046  0.014441  0.701326    0.259400
    # E12  GALILEO-1  0.366780  0.310270  0.195602  0.039986  0.945892    0.765318
    # E19  GALILEO-1  0.154111  0.141690  0.060615  0.013444  0.284842    0.244182
    for field in fields:
        extra_row_names = list()
        extra_rows = list()
        df_field = pd.DataFrame(columns=columns)

        # Determine functional values for each satellite
        for satellite in sorted(dset.unique("satellite")):
            row = [_get_satellite_type(dset, satellite)]
            for function in functions:
                try:
                    value = dset.apply(function, field, satellite=satellite)
                except:
                    value = float("nan")
                row.append(value)
            df_row = pd.DataFrame([row], columns=columns, index=[satellite])
            df_field = df_field.append(df_row)

        # Sort dataframe after satellite type -> TODO: Better solution for sorting after index?
        df_field["satellite"] = df_field.index
        df_field = df_field.sort_values(by=["type", "satellite"])
        del df_field["satellite"]

        # Determine functional values for each system
        for system in sorted(dset.unique("system")):
            row = [""]  # Append satellite type
            for function in functions:
                try:
                    value = dset.apply(function, field, system=system)
                except:
                    value = float("nan")
                row.append(value)
            system_name = f"__SYSTEM_{system}__"
            extra_row_names.append(system_name)
            extra_rows.append(pd.DataFrame([row], columns=columns, index=[system_name]))

        # Determine functional values for each satellite type
        for type_ in sorted(dset.unique("satellite_type")):
            row = [""]  # Append satellite type
            for function in functions:
                try:
                    value = dset.apply(function, field, satellite_type=type_)
                except:
                    value = float("nan")
                row.append(value)
            type_name = f"__{type_}__"
            extra_row_names.append(type_name)
            extra_rows.append(pd.DataFrame([row], columns=columns, index=[type_name]))

        # Append extra rows
        for row in extra_rows:
            df_field = df_field.append(row)
        df_field = df_field.reindex(
            columns=columns
        )  # TODO: Why is the column order be changed by appending extra rows?

        # Add field Dataframe to field Dataframe dictionary
        field_dfs.update({field: df_field})

    _plot_bar_stacked_sisre(fid, field_dfs, extra_row_names, figure_dir)
    _plot_bar_stacked_sisre_satellites(fid, field_dfs, extra_row_names, figure_dir)

    # Write field Dataframes
    for field, df in field_dfs.items():
        column = "rms"  # TODO: Loop over columns?
        fid.write("\n\n#{} RMS for all satellites in [m]\n".format(field.upper()))
        _write_dataframe_to_markdown(fid, df, float_format="6.3f")
        _plot_bar_dataframe_columns(fid, figure_dir, df, field, extra_row_names, column=column)


def _statistics(fid, figure_dir, dset):
    _statistics_satellite(fid, figure_dir, dset)


def _unhealthy_satellites(fid, dset):
    """Write overview over unhealthy satellites

    Args:
       fid:
       dset:
    """
    brdc = apriori.get(
        "orbit",
        rundate=dset.rundate,
        time=dset.time,
        satellite=tuple(dset.satellite),
        system=tuple(dset.system),
        station=dset.dataset_name.upper(),
        apriori_orbit="broadcast",
    )
    fid.write("{:25s} = {sats:60s}\n" "".format("Unhealthy satellites", sats=" ".join(brdc.unhealthy_satellites())))


def _write_config(fid):
    """Print the configuration options for a given technique and model run date

    Args:
       fid:
    """
    # Print the individual configuration options
    fid.write("#Configuration of SISRE analysis\n__Located at {}__\n".format(", ".join(config.tech.sources)))
    fid.write("```\n")
    fid.write(str(config.tech.as_str(key_width=25, width=70, only_used=True)))
    fid.write("\n```\n")
    fid.write("\\newpage\n")


def _write_dataframe_to_markdown(fid, df, float_format=""):
    """Write Pandas DataFrame to Markdown

    Args:
        fid (_io.TextIOWrapper):    File object
        df (DataFrame):             Pandas DataFrame
        float_format (str):         Define formatters for float columns
    """
    column_length = [len(c) for c in df.columns]

    # Write header
    if list(df.index):  # Add DataFrame index to header
        num_space = len(max(df.index))
        head_line_1 = "\n| {} ".format(" " * num_space)
        head_line_2 = "|-{}-".format("-" * num_space)
    else:
        header_1 = ""

    fid.write(head_line_1 + "| {} |\n".format(" | ".join(list(df.columns))))
    fid.write(head_line_2 + "|-{}:|\n".format("-|-".join([n * "-" for n in column_length])))

    # Write data
    for index, row in df.iterrows():

        line = "| {idx:s} |".format(idx=index) if index else ""  # Add DataFrame index column

        for _, v in row.items():
            if isinstance(v, float):
                line = line + " {:{fmt}} |".format(v, fmt=float_format)
            else:
                line = line + " {} |".format(v)
        fid.write(line + "\n")
    fid.write("\n")
