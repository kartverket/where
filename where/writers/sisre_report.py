"""Write report about a SISRE analysis run (only for one session and not several sessions (stations))

Description:
------------

asdf




$Revision: 15265 $
$Date: 2018-06-05 23:34:59 +0200 (Tue, 05 Jun 2018) $
$LastChangedBy: dahmic $
"""
# Standard library imports
from collections import namedtuple
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
from where.lib import files
from where.lib import gnss
from where.lib import plugins
from where.reports import report


SubplotConfig = namedtuple("SubplotConfig", ["ylabel", "color", "ydata"])
SubplotConfig.__doc__ = """A convenience class for defining a field for subplot configuration

    Args:
        ylabel (str):           Y-axis label
        color (str):            Color of scatter plot
        ydata (numpy.ndarray):  Y-axis data
    """

@plugins.register
def sisre_report(dset):
    """Write SISRE report

    Args:
        dset:       Dataset, a dataset containing the data.
    """
    with files.open("output_gnss_session_report", mode="wt") as fid:
        _header(fid)
        _write_config(fid)
        fid.write("\n# Satellite status\n\n")
        _unhealthy_satellites(fid, dset)
        _eclipse_satellites(fid, dset)

        # Generate figure directory to save figures generated for SISRE report
        figure_dir = files.path("output_sisre_report_figure", file_vars=config.files.vars)
        figure_dir.mkdir(parents=True, exist_ok=True)
        _plot_satellite_bias(fid, figure_dir, dset)
        _plot_orbit_and_clock_differences(fid, figure_dir, dset)
        _plot_sisre(fid, figure_dir, dset)
        _statistics(fid, figure_dir, dset)


def _eclipse_satellites(fid, dset):
    """Write overview over satellites in eclipse

    Args:
       fid:
       dset:
    """
    #TODO: time period of eclipting satellites needed
    brdc = apriori.get(
        "orbit",
        rundate=dset.rundate,
        time=dset.time,
        satellite=tuple(dset.satellite),
        system=tuple(dset.system),
        station=dset.dataset_name.upper(),
        apriori_orbit="broadcast",
    )

    fid.write("{:25s} = {sats:47s}\n" "".format('Satellites in eclipse', sats=" ".join(gnss.check_satellite_eclipse(brdc.dset))))


def _header(fid):
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


def _get_satellite_type(dset, satellite):
        idx = dset.filter(satellite=satellite)
        return set(dset.satellite_type[idx]).pop()

def _plot_bar_dataframe_columns(fid, figure_dir, df, field, extra_row_names, column='rms'):
    """Generate bar plot of given dataframe columns (colored and ordered by satellite type)
    """
    fontsize=12
    df_reduced = df.drop(extra_row_names)  # Remove extra rows
    
    # Assign to each satellite type a color
    colors = dict() 
    #TODO: Better handling of color definition?   
    #color_def = ['cornflowerblue', 'firebrick', 'violet', 'gold', 'limegreen', 'deepskyblue', 'orangered']
    color_def = ['red', 'tomato', 'lightsalmon', 'blue', 'royalblue', 'deepskyblue', 'paleturquoise']
    #color_def = ['C'+str(idx) for idx in range(0,10)]
    for type_ in sorted(set(df_reduced.type)):
        colors.update({type_: color_def.pop()})

    # Generate bar plot
    df_color = df_reduced['type'].apply(lambda x: colors[x])
    fig_width = len(df_reduced.index)/4 if len(df_reduced.index) > 30 else 6.4
    ax = df_reduced[column].plot(kind="bar", color=df_color, width=0.8, figsize=(fig_width, fig_width/1.33))
    ax.set_xlabel("Satellite", fontsize=fontsize)
    ax.set_ylabel(f'{field.upper()} {column.upper()} [m]', fontsize=fontsize)

    # Make legend
    satellite_type_patch = [mpatches.Patch(color=v, label=k) for k, v in sorted(colors.items())]
    ax.legend(handles=satellite_type_patch, bbox_to_anchor=(1.04, 1), loc=2, borderaxespad=0., ncol=1)

    plt.tight_layout()
    plt.savefig(figure_dir/f'plot_bar_{field}_{column}.pdf')
    plt.clf() # clear the current figure

    fid.write(f'![{field.upper()} {column.upper()} for all satellites sorted by satellite type]({figure_dir}/plot_bar_{field}_{column}.pdf)\n')
    fid.write('\\newpage\n')



def _plot_bar_stacked_sisre(fid, field_dfs, extra_row_names, figure_dir):
    figure_path = figure_dir/f'plot_bar_stacked_sisre_extra_rows.pdf'

    # Remove satellite rows
    df_sisre_extra = field_dfs['sisre'].drop(set(field_dfs['sisre'].index) - set(extra_row_names))
    df_sisre_orb_extra = field_dfs['sisre_orb'].drop(set(field_dfs['sisre_orb'].index) - set(extra_row_names))

    _plot_bar_stacked(fid, df_sisre_extra, df_sisre_orb_extra, figure_path=figure_path, xticks_rotation=20)

    fid.write(f'![Blue bars indicate orbit-only SISRE RMS, orange bars clock-only SISRE RMS and green bars 95th percentile SISRE]({figure_path})\n')
    fid.write('\\newpage\n')


def _plot_bar_stacked_sisre_satellites(fid, field_dfs, extra_row_names, figure_dir):
    figure_path = figure_dir/f'plot_bar_stacked_sisre.pdf'

    # Remove extra rows
    df_sisre = field_dfs['sisre'].drop(extra_row_names)
    df_sisre_orb = field_dfs['sisre_orb'].drop(extra_row_names)

    _plot_bar_stacked(fid, df_sisre, df_sisre_orb, figure_path=figure_path, xlabel="Satellite")

    fid.write(f'![Blue bars indicate orbit-only SISRE RMS, orange bars clock-only SISRE RMS and green bars 95th percentile SISRE]({figure_path})\n')
    fid.write('\\newpage\n')


def _plot_bar_stacked(fid, df_sisre, df_sisre_orb, figure_path, xlabel='', xticks_rotation=None):
    """Generate bar plot of given dataframe columns (colored and ordered by satellite type)
    """
    fontsize=12

    # Generate new dataframe with columns 'orbit-only SISRE' and 'clock-only SISRE'
    data = np.hstack((np.array([df_sisre_orb.rms]).T, (np.array([df_sisre.rms]) - np.array([df_sisre_orb.rms])).T, (np.array([df_sisre.percentile]) - np.array([df_sisre.rms])).T)) 
    df_merged = pd.DataFrame(data=data, index=df_sisre.index, columns=['orbit-only SISRE RMS', 'clock-only SISRE RMS', '95th percentile SISRE'])

    # Generate bar plot
    fig_width = len(df_merged.index)/4 if len(df_merged.index) > 30 else 6.4
    ax = df_merged.plot(kind="bar", stacked=True, width=0.8, figsize=(fig_width, fig_width/1.33))
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(f'SISRE [m]', fontsize=fontsize)
    #ax.legend(bbox_to_anchor=(1.04, 1), loc=2, borderaxespad=0., ncol=1)
    #plt.axhline(y=2, linewidth=2, color='#d62728') #TODO: Dependet on GNSS!!!
    ax.legend(bbox_to_anchor=(1, 1.15), loc=1, borderaxespad=0., ncol=2)

    if xticks_rotation is not None:
        plt.xticks(rotation=xticks_rotation)
    plt.tight_layout()
    plt.savefig(figure_path)
    plt.clf() # clear the current figure



def _plot_satellite_bias(fid, figure_dir, dset):

    for field, orbit in {'bias_brdc': 'Broadcast', 'bias_precise': 'Precise'}.items():
        if np.sum(dset[field]) != 0:
            plt.scatter(dset.time.gps.datetime,dset[field])
            plt.ylabel(f'{orbit} satellite bias [m]')
            plt.xlim([min(dset.time.gps.datetime),  max(dset.time.gps.datetime)])
            plt.xlabel('Time [GPS]')
            plt.savefig(figure_dir/f'plot_scatter_{field}.pdf')

            fid.write(f'![{orbit} satellite bias]({figure_dir}/plot_scatter_{field}.pdf)\n')
            fid.write('\\newpage\n')


def _plot_orbit_and_clock_differences(fid, figure_dir, dset):
    
    # Define configuration of subplots
    subplots = (
    SubplotConfig("Δalong-track [m]", "paleturquoise", dset.orb_diff.itrs[:, 0]),
    SubplotConfig("Δcross-track [m]", "deepskyblue", dset.orb_diff.itrs[:, 1]),
    SubplotConfig("Δradial [m]", "royalblue", dset.orb_diff.itrs[:, 2]),
    SubplotConfig("Δclock [m]", "tomato", dset.clk_diff_sys),
    )

    _plot_scatter_subplots(dset.time.gps.datetime, subplots, figure_path=figure_dir/f'plot_scatter_orbit_clock_differences.pdf', xlabel='Time [GPS]' )   

    fid.write(f'![Broadcast ephemeris errors for along-track, cross-track and radial orbit errors and clock errors]({figure_dir}/plot_scatter_orbit_clock_differences.pdf)\n')
    fid.write('\\newpage\n')


def _plot_sisre(fid, figure_dir, dset):
    
    # Define configuration of subplots
    subplots = (
    SubplotConfig("orbit-only SISRE [m]", "paleturquoise", dset.sisre_orb),
    SubplotConfig("clock-only SISRE [m]", "deepskyblue", dset.sisre - dset.sisre_orb),
    SubplotConfig("SISRE [m]", "royalblue", dset.sisre),
    )

    _plot_scatter_subplots(dset.time.gps.datetime, subplots, figure_path=figure_dir/f'plot_scatter_sisre.pdf', xlabel='Time [GPS]' )   

    fid.write(f'![Orbit-only SISRE, clock-only SISRE and SISRE]({figure_dir}/plot_scatter_sisre.pdf)\n')
    fid.write('\\newpage\n')


def _plot_scatter_subplots(xdata, subplots, figure_path, xlabel=''):
    marker = "." #point marker type

    fig, axes = plt.subplots(len(subplots),1, sharex=True, sharey=True)
    fig.set_figheight(9) #inches
    for idx, ax in enumerate(axes):
        ax.set(ylabel=subplots[idx].ylabel)
        ax.set_xlim([min(xdata), max(xdata)]) # otherwise time scale of x-axis is not correct -> Why?
        text = f'mean $ = {np.mean(subplots[idx].ydata):.2f} \pm {np.std(subplots[idx].ydata):.2f}$ m'
        ax.text(0.98, 0.98, text, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)                      
        ax.scatter(xdata, subplots[idx].ydata, marker=marker, color=subplots[idx].color)
        ax.set(xlabel=xlabel)

    fig.autofmt_xdate() # rotates and right aligns the x labels, and moves the bottom of the
                        # axes up to make room for them
    plt.tight_layout()
    plt.savefig(figure_path)
    plt.clf() # clear the current figure


    



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
    columns   = ["type", "rms", "mean", "std", "min", "max", "percentile"]
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
        for satellite in sorted(dset.unique('satellite')):   
            row = [_get_satellite_type(dset, satellite)]
            for function in functions:
                try:
                    value = dset.apply(function, field, satellite=satellite)
                except:
                    value = float('nan')
                row.append(value)
            df_row = pd.DataFrame([row], columns=columns, index=[satellite])
            df_field = df_field.append(df_row)

        # Sort dataframe after satellite type -> TODO: Better solution for sorting after index?
        df_field["satellite"] = df_field.index
        df_field = df_field.sort_values(by=["type", "satellite"])
        del df_field["satellite"]


        # Determine functional values for each system
        for system in sorted(dset.unique('system')):
            row = ['']  # Append satellite type
            for function in functions:
                try:
                    value = dset.apply(function, field, system=system)
                except:
                    value = float('nan')
                row.append(value)
            system_name = f'__SYSTEM_{system}__'
            extra_row_names.append(system_name)
            extra_rows.append(pd.DataFrame([row], columns=columns, index=[system_name]))
            
        # Determine functional values for each satellite type
        for type_ in sorted(dset.unique('satellite_type')):
            row = ['']  # Append satellite type
            for function in functions:
                try:
                    value = dset.apply(function, field, satellite_type=type_)
                except:
                    value = float('nan')
                row.append(value)
            type_name = f'__{type_}__'
            extra_row_names.append(type_name)
            extra_rows.append(pd.DataFrame([row], columns=columns, index=[type_name]))

        # Append extra rows 
        for row in extra_rows:
            df_field = df_field.append(row)
        df_field = df_field.reindex(columns=columns) #TODO: Why is the column order be changed by appending extra rows?

        # Add field Dataframe to field Dataframe dictionary
        field_dfs.update({field: df_field})

    _plot_bar_stacked_sisre(fid, field_dfs, extra_row_names, figure_dir)
    _plot_bar_stacked_sisre_satellites(fid, field_dfs, extra_row_names, figure_dir)

    # Write field Dataframes
    for field, df in field_dfs.items():
        column = 'rms' #TODO: Loop over columns?
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
    fid.write("{:25s} = {sats:60s}\n" "".format('Unhealthy satellites', sats=" ".join(brdc.unhealthy_satellites())))

def _write_config(fid):
    """Print the configuration options for a given technique and model run date

    Args:
       fid:
    """
    # Print the individual configuration options
    fid.write("#Configuration of SISRE analysis\n__Located at {}__\n".format(", ".join(config.tech.sources)))
    fid.write("```\n")
    fid.write(
        str(config.tech.as_str(key_width=25, width=70, only_used=True))
    )
    fid.write("\n```\n")

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
