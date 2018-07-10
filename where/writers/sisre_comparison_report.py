"""Compare different SISRE Where datasets

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
from where.lib import log
from where.lib import gnss
from where.lib import plugins


@plugins.register
def sisre_comparison_report(dset):
    """Compare SISRE datasets

    Args:
        dset (list):       List with different SISRE datasets. The datasets contain the data.
    """
    dsets = dset
    df_merged = pd.DataFrame()

    for name, dset in dsets.items():
        user_type_name = _get_user_type_name(name)
        df = dset.as_dataframe(fields=["satellite", "system", "sisre", "time.gps"])  # , index="time.gps")
        df = df.rename(columns={"sisre": user_type_name})
        if df_merged.empty:
            df_merged = df
            continue
        df_merged = df_merged.merge(df, on=["satellite", "system", "time.gps"], how="outer")

    with files.open(
        file_key="output_sisre_comparison_report", file_vars=dsets[next(iter(dsets))].vars, mode="wt"
    ) as fid:
        # _header(fid)
        fid.write("#SISRE analysis\n")

        # Generate figure directory to save figures generated for SISRE report
        figure_dir = files.path("output_sisre_report_figure", file_vars=dset.vars)
        figure_dir.mkdir(parents=True, exist_ok=True)

        # _plot_bar_sisre_satellite_percentile(df_merged, fid, figure_dir)
        _plot_bar_sisre_constellation_percentile(df_merged, fid, figure_dir)
        exit()  # MURKS

        # if write_level <= enums.get_value("write_level", "detail"):
        #    _plot_scatter_satellite_bias(fid, figure_dir, dset)
        #    _plot_scatter_orbit_and_clock_differences(fid, figure_dir, dset)
        # _plot_histogram_sisre(fid, figure_dir, dset)
        # _plot_scatter_sisre(fid, figure_dir, dset)
        # _statistics(fid, figure_dir, dset)


def _get_user_type_name(dset_name):
    # TODO: Pattern search inav_e1 does not work. How to solve?
    user_type_name_def = {"fnav_e1e5a": "E1E5a", "inav_e1_": "E1", "inav_e1e5b": "E1E5b", "lnav_l1l2": "L1L2"}
    user_type_name = None

    for pattern, type_name in user_type_name_def.items():
        if pattern in dset_name:
            return type_name

    if not user_type_name:
        log.fatal(f"User type name not found in {dset_name}.")


def _plot_bar_sisre_constellation_percentile(df, fid, figure_dir):
    # df_monthly_percentile = df.set_index('time.gps').resample('D', how=lambda x: np.nanpercentile(x, q=95))
    df_monthly_percentile = df.dropna().set_index("time.gps").resample("M", how=lambda x: np.percentile(x, q=95))
    df_monthly_percentile.index = df_monthly_percentile.index.strftime("%b-%Y")
    df_monthly_percentile.transpose().plot(kind="bar")
    plt.axhline(2, color="r")
    plt.xlabel("Signal combination for single- and dual-frequency users")
    plt.xticks(rotation=0)
    plt.ylabel("SISRE [m] (95th percentile)")
    plt.legend(bbox_to_anchor=(1.04, 1), loc=2, borderaxespad=0., ncol=1)
    # plt.legend(bbox_to_anchor=(1, 1.15), loc=1, borderaxespad=0., ncol=3)
    plt.tight_layout()
    plt.savefig(figure_dir / f"plot_bar_sisre_constellation_percentile.pdf")
    plt.clf()  # clear the current figure

    fid.write(
        f"![Monthly 95th percentile of global average SISRE for single- and dual-frequency users]({figure_dir}/plot_bar_sisre_constellation_percentile.pdf)\n"
    )
    fid.write("\\newpage\n")


def _plot_bar_sisre_satellite_percentile(df, fid, figure_dir):
    import IPython

    IPython.embed()

    user_types = ["E1", "E1E5a", "E1E5b"]

    fig, axes = plt.subplots(len(user_types), 1, sharex=True, sharey=True)  # figsize=(6, 6));
    # plt.subplots_adjust(wspace=0.5, hspace=0.5);
    fig.set_figheight(9)  # inches
    # fig.suptitle(f"{title}", y=1.0)
    for idx, (ax, user) in enumerate(zip(axes, user_types)):

        df_user = df.pivot(index="time.gps", columns="satellite", values=user)
        df_user_monthly_percentile = df_user.dropna().resample("M", how=lambda x: np.percentile(x, q=95))
        df_user_monthly_percentile.index = df_user_monthly_percentile.index.strftime("%b-%Y")
        # df_user_monthly_percentile.transpose().plot(kind="bar")
        # df_user_monthly_percentile.transpose().plot(subplots=True, kind="bar", ax=ax, legend=False, sharex=False, sharey=False);
        plt.bar

        # df_merged.plot(kind="bar")

        # ax.set(ylabel="95th percentile SISRE [m]")
        # ax.set_xlim([min(xdata), max(xdata)])  # otherwise time scale of x-axis is not correct -> Why?
        # text = f"mean $ = {np.mean(subplots[idx].ydata):.2f} \pm {np.std(subplots[idx].ydata):.2f}$ m"
        # ax.text(0.98, 0.98, text, horizontalalignment="right", verticalalignment="top", transform=ax.transAxes)
        # ax.bar()
        # ax.set(xlabel='Satellite')

    for user in user_types:
        df_user = df.pivot(index="time.gps", columns="satellite", values=user)
        df_user_monthly_percentile = df_user.dropna().resample("M", how=lambda x: np.percentile(x, q=95))
        df_user_monthly_percentile.index = df_user_monthly_percentile.index.strftime("%b-%Y")
        df_user_monthly_percentile.transpose().plot(kind="bar")
        plt.axhline(7, color="r")
        plt.xlabel("Satellite")
        plt.ylabel("SISRE [m]")
        # plt.legend(bbox_to_anchor=(1.04, 1), loc=2, borderaxespad=0., ncol=1)
        plt.legend(bbox_to_anchor=(1, 1.15), loc=1, borderaxespad=0., ncol=3)

        plt.show()

    plt.tight_layout()
    plt.savefig(figure_dir / f"plot_bar_sisre_constellation_percentile.pdf")
    plt.clf()  # clear the current figure

    fid.write(
        f"![Monthly 95th percentile of global average SISRE for single- and dual-frequency users]({figure_dir}/plot_bar_sisre_constellation_percentile.pdf)\n"
    )
    fid.write("\\newpage\n")


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
  \\label{eq:sisre}  df
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
def _plot_histogram_sisre(fid, figure_dir, dset):
    # import IPython; IPython.embed()
    for sys in dset.unique("system"):
        idx = dset.filter(system=sys)
        mean = np.mean(dset.sisre[idx])
        std = np.std(dset.sisre[idx])
        plt.hist(dset.sisre[idx], normed=True, bins=30)
        plt.text(0.98, 0.98, f"$mean={mean:5.3f}\ \pm {std:5.3f}$ m")
        # text(0.98, 0.98, text, horizontalalignment="right", verticalalignment="top", transform=ax.transAxes)
        plt.xlabel("SISRE [m]")
        plt.ylabel("Frequency")
        # plt.show()


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
