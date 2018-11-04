"""Write report about a SISRE analysis run, whereby several stations with different navigation files can be processed.

Description:
------------

Normally the SISRE analysis is based only on one session/station. But in principle the SISRE analysis can be used for several station navigation files. In this case the session option in the configuration file is used for defining a list of stations, which are processed within the SISRE analysis. 

"""
# Standard library imports
from datetime import datetime
import re
import textwrap

# External library imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# WHERE imports
import where
from where import apriori
from where import cleaners
from where.lib import config
from where.lib import files
from where.lib import gnss
from where.lib import plugins
from where.reports import report


@plugins.register
def write_session_report(rundate, tech):
    dsets = dict()
    for _, dset, report_data in report.reports("remover_data"):
        remover_name = report_data["remover_name"]
        station = dset.dataset_name
        dset_station = dsets.setdefault("removers", dict()).setdefault(station, dset)
        dset_station.add_float(
            "keep_idx_{}".format(remover_name), val=report_data["keep_idx"]
        )  # TODO: What is keep_idx?

    for _, dset, report_data in report.reports("orbit_data"):
        dsets.setdefault("orbit", dict())[report_data["station"]] = dset

    with files.open("output_gnss_session_report", mode="wt") as fid:
        header(fid)
        write_config(fid)
        unhealthy_satellites(fid, dsets)
        eclipse_satellites(fid, dsets)
        # rejected_satellites_per_station(fid, dsets)
        # rejected_satellite_observations(fid, dsets)
        # rejected_station_observations(fid, dsets)
        statistics(fid, dsets)


def header(fid):
    _, _, title = next(report.reports("title"))

    fid.write(
        """---
title: {title}
author: Where v{version} [{user}]
date: {nowdate:%Y-%m-%d}
---
""".format(
            title=title["text"], version=where.__version__, user=config.analysis.user.str, nowdate=datetime.now()
        )
    )

    fid.write("\\newpage\n")
    fid.write("#SISRE analysis\n")

    fid.write("\\begin{equation}\n")
    fid.write(
        "     \\text{SISRE(orb)}^s = \sqrt{w_r^2 \cdot \Delta {r^s}^2 + w_{a,c}^2 \cdot (\Delta {a^s}^2 + \Delta {c^s}^2)}\n"
    )
    fid.write("  \\label{eq:sisre_orb}\n")
    fid.write("\\end{equation}\n")
    fid.write("\n")
    fid.write("\\begin{equation}\n")
    fid.write(
        "     \\text{SISRE}^s = \sqrt{(w_r \cdot \Delta r^s - \Delta t^s)^2 + w_{a,c}^2 \cdot (\Delta {a^s}^2 + \Delta {c^s}^2)}\n"
    )
    fid.write("  \\label{eq:sisre}\n")
    fid.write("\\end{equation}\n")
    fid.write("\n")
    fid.write("\\begin{tabular}{lll}\n")
    fid.write("with\\\n")
    fid.write(" & $\\text{SISRE(orb)}^s$  & Orbit-only SISRE for satellite $s$, \\\\\n")
    fid.write(" & $\\text{SISRE}^s$       & SISRE for satellite $s$, \\\\\n")
    fid.write(
        " &$\Delta a^s$, $\Delta c^s$, $\Delta r^s$   & Satellite coordinate differences between broadcast and precise \\\\\n"
    )
    fid.write(" &                        & ephemeris in along-track, cross-track and radial for satellite $s$,\\\\\n")
    fid.write(
        " &$w_r$                   & SISRE weight factor for radial errors (see Table \\ref{tab:sisre_weight_factors}),\\\\\n"
    )
    fid.write(" &$w_{a,c}$               & SISRE weight factor for along-track and cross-track errors \\\\\n")
    fid.write(" &                        & (see Table \\ref{tab:sisre_weight_factors}).\\\\\n")
    fid.write("\end{tabular}\n")
    fid.write("\n")
    fid.write("\\begin{table}[!ht]\n")
    fid.write("\\begin{center}\n")
    fid.write("  \\begin{tabular}[c]{lll}\n")
    fid.write("    \hline\n")
    fid.write("    System       & $w_r$  & $w_{a,c}^2$ \\\\\n")
    fid.write("    \hline\n")
    fid.write("    GPS          & $0.98$ & $1/49$ \\\\\n")
    fid.write("    GLONASS      & $0.98$ & $1/45$ \\\\\n")
    fid.write("    Galileo      & $0.98$ & $1/61$ \\\\\n")
    fid.write("    QZSS         & $0.99$ & $1/126$ \\\\\n")
    fid.write("    BeiDou (MEO) & $0.98$ & $1/54$ \\\\\n")
    fid.write("    \hline\n")
    fid.write("  \end{tabular}\n")
    fid.write(
        "  \caption[SISRE weight factors for radial ($w_r$) and along-track and cross-track errors ($w_{a,c}$)]{SISRE weight factors for radial ($w_r$) and along-track and cross-track ($w_{a,c}$) errors }\n"
    )
    fid.write("  \label{tab:sisre_weight_factors}\n")
    fid.write("\end{center}\n")
    fid.write("\end{table}\n")
    fid.write("\\newpage\n")


def write_config(fid):
    """Print the configuration options for a given technique and model run date

    Args:
       fid:
    """
    skip_config_keys = [
        "reference_ellipsoid",
        "reference_frame",
        "ocean_tides",
        "atmospheric_tides",
        "mean_pole_version",
        "eop_models",
    ]  # TODO

    # Print the individual configuration options
    fid.write("#Configuration of SISRE\n__Located at {}__\n".format(", ".join(config.tech.sources)))
    fid.write("```\n")
    # TODO hjegei: Only_used option does not show all used options. Why?
    fid.write(
        str(config.tech.as_str(key_width=25, width=70, only_used=True))
    )  # TODO: extra str() should not be necessary?
    fid.write("\n```\n")


def unhealthy_satellites(fid, dsets):
    """Write overview over unhealthy satellites

    Args:
       fid:
       dsets:
    """
    stations = sorted(dsets["orbit"].keys())
    if not stations:
        return 0

    fid.write("\n# Overview over unhealthy satellites\n\n")
    fid.write("| Sta.  |   Satellites                             |\n")
    fid.write("|-------|-----------------------------------------:|\n")

    for station in stations:
        dset = dsets["orbit"][station]

        brdc = apriori.get(
            "orbit",
            rundate=dset.rundate,
            time=dset.time,
            satellite=tuple(dset.satellite),
            system=tuple(dset.system),
            station=dset.dataset_name.upper(),
            apriori_orbit="broadcast",
        )
        fid.write("| {sta:5s} | {sats:40s} |\n" "".format(sta=station, sats=" ".join(brdc.unhealthy_satellites())))


def eclipse_satellites(fid, dsets):
    """Write overview over satellites in eclipse

    Args:
       fid:
       dsets:
    """
    stations = sorted(dsets["orbit"].keys())
    if not stations:
        return 0

    fid.write("\n# Overview over satellites in eclipse\n\n")
    fid.write("| Satellites                                      |\n")
    # fid.write('| Sat. | From    | To                             |\n') #TODO: time period of eclipting satellites needed
    fid.write("|------------------------------------------------:|\n")

    dset = dsets["orbit"][stations[0]]
    brdc = apriori.get(
        "orbit",
        rundate=dset.rundate,
        time=dset.time,
        satellite=tuple(dset.satellite),
        system=tuple(dset.system),
        station=dset.dataset_name.upper(),
        apriori_orbit="broadcast",
    )
    fid.write("| {sats:47s} |\n" "".format(sats=" ".join(gnss.check_satellite_eclipse(brdc.dset))))


# def configuration(fid, dsets):
#
#    cfg = config.get_tech_config()
#    config_list = ['sampling_rate', 'systems']
#
#    fid.write('# Configuration settings\n')
#    for option in config_list:
#        fid.write('{option:14s} = {value:20s}\n'.format(option=option, value=config.get(option, cfg)))
#
#    #satellites = config.getlist('ignore_satellite', cleaners.get_config(), section=dset.dataset_name)


def rejected_satellites_per_station(fid, dsets):

    for dset in dsets.values():
        fid.write(
            "\n# Overview over rejected observations for all satellites for station {}\n\n"
            "".format(dset.dataset_name.upper())
        )
        fid.write("| Sat.  |   #obs |  #eobs | #eobs% |\n")
        fid.write("|-------|-------:|-------:|-------:|\n")
        for sat in dset.unique("satellite"):
            sat_idx = dset.filter(satellite=sat)
            num_obs = sum(sat_idx)
            num_edit_obs = sum(np.logical_not(dset.keep_idx[sat_idx]))
            fid.write(
                "| {sat:5s} | {obs:6d} | {eobs:6d} | {eobsp:6.0f} |\n"
                "".format(sat=sat, obs=num_obs, eobs=num_edit_obs, eobsp=100 * num_edit_obs / num_obs)
            )

        num_edit_obs = sum(np.logical_not(dset.keep_idx))
        fid.write(
            "| Sum   | {obs:6d} | {eobs:6d} | {eobsp:6.0f} |\n"
            "".format(obs=dset.num_obs, eobs=num_edit_obs, eobsp=100 * num_edit_obs / dset.num_obs)
        )


def rejected_satellite_observations(fid, dsets):

    sum_num_obs = sum_num_edit_obs = 0
    fid.write("\n# Overview over rejected observations for all satellites\n\n")
    fid.write("|   Sat.    |     #obs   |    #eobs   |   #eobs%   |\n")
    fid.write("|-----------|-----------:|-----------:|-----------:|\n")
    for dset in dsets.values():
        sum_num_obs += dset.num_obs
        num_edit_obs = sum(np.logical_not(dset.keep_idx))
        sum_num_edit_obs += num_edit_obs
        fid.write(
            "|   {sta:5s}   |   {obs:6d}   |   {eobs:6d}   |   {eobsp:6.0f}   |\n"
            "".format(
                sta=dset.dataset_name, obs=dset.num_obs, eobs=num_edit_obs, eobsp=100 * num_edit_obs / dset.num_obs
            )
        )
    fid.write(
        "| __Sum__   | __{obs:6d}__ | __{eobs:6d}__ | __{eobsp:6.0f}__ |\n"
        "".format(obs=sum_num_obs, eobs=sum_num_edit_obs, eobsp=100 * sum_num_edit_obs / sum_num_obs)
    )


def rejected_station_observations(fid, dsets):
    stations = sorted(dsets.keys())
    if not stations:
        return

    satellites = dsets[stations[0]].unique("satellite")  # TODO: Better?
    fid.write("\n# Overview over rejected observations for all stations\n\n")
    fid.write("| Sta.  |   #obs |  #eobs | #eobs% |\n")
    fid.write("|-------|-------:|-------:|-------:|\n")
    for sat in satellites:
        sum_num_obs = sum_num_edit_obs = 0
        for station in stations:
            dset = dsets[station]
            sat_idx = dset.filter(satellite=sat)
            num_obs = sum(sat_idx)
            sum_num_obs += num_obs
            num_edit_obs = sum(np.logical_not(dset.keep_idx[sat_idx]))
            sum_num_edit_obs += num_edit_obs
        fid.write(
            "| {sat:5s} | {obs:6d} | {eobs:6d} | {eobsp:6.0f} |\n"
            "".format(sat=sat, obs=sum_num_obs, eobs=sum_num_edit_obs, eobsp=100 * sum_num_edit_obs / sum_num_obs)
        )


def statistics(fid, dsets):
    # TODO: sort satellites after satellite type
    _statistics_station_satellite(fid, dsets)


def _statistics_station_satellite(fid, dsets):
    """Generate statistics for each field, station and satellite

    Args:
        fid (_io.TextIOWrapper):    File object
        dsets (list):               List of Datasets
    """
    fields = ["sisre", "sisre_orb", "orb_diff_3d", "clk_diff"]
    stations = sorted(list(dsets["orbit"].keys()))
    columns = ["system"] + stations
    field_dfs = dict()
    satellites = set()

    # Get a unique set of satellites, which are valid for all stations
    for station in dsets["orbit"]:
        satellites = satellites.union(set(dsets["orbit"][station].unique("satellite")))

    # Generate field DataFrames with determined RMS values and with the satellites as indices and stations as columns
    #
    # Example:
    #      stas   trds
    # G01  1.589  1.465
    # G02  2.489  2.450
    # G03  1.234  1.197
    for field in fields:
        extra_rows = list()
        df_field = pd.DataFrame(columns=columns)

        # Determine RMS value for each satellite for given field and station
        for satellite in sorted(satellites):
            row = [satellite[0]]  # Get system identifier

            for station in sorted(stations):
                dset = dsets["orbit"][station]
                # TODO: What kind of exception?
                try:
                    rms = dset.rms(field, satellite=satellite)
                except:
                    rms = float(nan)
                row = row + [rms]
            df_row = pd.DataFrame([row], columns=columns, index=[satellite])
            df_field = df_field.append(df_row)

        # Add a MEAN column to the DataFrame
        df_field["MEAN"] = df_field.mean(axis=1)

        # Determine mean values for GNSS specific satellites and for all satellites
        # TODO: Determine RMS instead of MEAN
        for sys in df_field["system"].unique():
            mean_sys = df_field[(df_field.system == sys)].mean(axis=0)
            mean_sys.name = "__MEAN_{}__".format(sys)
            extra_rows.append(mean_sys)

        mean_row = df_field.mean(axis=0)
        mean_row.name = "__MEAN__"
        extra_rows.append(mean_row)

        for row in extra_rows:
            df_field = df_field.append(row)

        # Add field Dataframe to field Dataframe dictionary
        field_dfs.update({field: df_field})

    # Write field Dataframes
    for field, df in field_dfs.items():
        fid.write("\n\n#{} RMS for all stations and satellites in [m]\n".format(field.upper()))
        _write_dataframe_to_markdown(fid, df, remove_columns=["system"], float_format="6.3f")
        # TODO: _plot_dataframe(fid, df, field, column='MEAN')


def _plot_dataframe(fid, df, field, column):
    # TODO: write it to file and add reference to report
    ax = df[column].plot(kind="bar")
    ax.set_xlabel("Satellite", fontsize=12)
    ax.set_ylabel(field + " [m]", fontsize=12)
    plt.show()


# TODO: write a routine, which calculates the RMS for a Pandas Dataframe (look after df.apply())
# def _df_rms(df, axis=0):


def _write_dataframe_to_markdown(fid, df, remove_columns=[], float_format=""):
    """Write Pandas DataFrame to Markdown

    Args:
        fid (_io.TextIOWrapper):    File object
        df (DataFrame):             Pandas DataFrame
        remove_columns (list):      List with columns to remove from DataFrame
        float_format (str):         Define formatters for float columns
    """
    if remove_columns:
        for column in remove_columns:
            del df[column]

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
