"""Compare different SISRE Where datasets

Description:
------------

A dictionary with datasets is used as input for this writer. The keys of the dictionary have to include the signal type
used for generation of the SISRE dataset. -> TODO: signal type of SISRE dataset has to be defined via meta variable!!!!


Example:
--------

    from where import data
    from where import writers

    # Read a dataset
    dset = data.Dataset(rundate=rundate, tech=tech, stage=stage, dataset_name=name, dataset_id=dataset_id)

    # Write dataset
    writers.write_one('sisre_comparison_report', dset=dsets, do_report=False)

"""
# Standard library imports
from pathlib import PosixPath
from typing import Any, Dict, List, Tuple, Union

# External library imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Midgard imports
from midgard.collections import enums
from midgard.dev import plugins

# Where imports
from where.lib import config, log
from where.writers._report import Report

_SECTION = "_".join(__name__.split(".")[-1:])

FIGURE_FORMAT = "png"
FILE_NAME = __name__.split(".")[-1]


@plugins.register
def sisre_comparison_report(dset: Dict[str, "Dataset"]) -> None:
    """Compare SISRE datasets
      
    Args:
        dset: Dictionary with SISRE solution name as keys (e.g. cnes_inav_e1, cnes_inav_e1e5b, cnes_fnav_e1e5a) and
              the belonging Dataset as value
    """
    dset_first = dset[list(dset.keys())[0]]
    file_vars = {**dset_first.vars, **dset_first.analysis}

    # Get file variables
    file_vars = {**dset_first.vars, **dset_first.analysis}
    for sol in dset.keys(): # Add service information (OS and/or HAS) to file variables
        if "HAS" in dset[sol].meta["service"]:
            file_vars["service"] = "has"
            file_vars["SERVICE"] = "HAS"
            break
        else:
            file_vars["service"] = "os"
            file_vars["SERVICE"] = "OS"

    # Generate figure directory to save figures generated for SISRE report
    figure_dir = config.files.path("output_sisre_comparison_report_figure", file_vars=file_vars)
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    df, df_month_perc, df_month_perc_rms, df_month_rms = _generate_dataframe(dset, field="sisre")
    _plot_bar_sisre_signal_combination_percentile(df_month_perc, figure_dir, threshold=False)
    _plot_bar_sisre_signal_combination_percentile(df_month_perc, figure_dir, threshold=True)
    _plot_bar_sisre_signal_combination_rms(df_month_rms, figure_dir)
    _plot_bar_sisre_signal_combination_percentile(df_month_perc_rms, figure_dir, threshold=False, suffix="_rms")
    _plot_bar_sisre_signal_combination_percentile(df_month_perc_rms, figure_dir, threshold=True, suffix="_rms")


    # Generate SISRE comparison report
    path = config.files.path("output_sisre_comparison_report", file_vars=file_vars)
    with config.files.open_path(path, create_dirs=True, mode="wt") as fid:

        rpt = Report(fid, rundate=dset_first.analysis["rundate"], path=path, description="Comparison of SISRE analyses")
        rpt.title_page()

        # +TODO: If possible should following lines added to _add_to_report().
        rpt.add_text("#Comparison of SISE analyses\n")
        rpt.add_text("In the following SISE analyses results are compared for:\n\n")
        rpt.add_text("* Monthly 95th percentile SISE for satellites\n")
        rpt.add_text("* Monthly 95th percentile and RMS SISE for signal combinations (users)\n")
        rpt.add_text("\\newpage\n")
        rpt.add_text("\n\n##Monthly 95th percentile SISE for satellites\n")
        # Produce plot with same yrange than for _plot_bar_sisre_signal_combination_percentile threshold plot
        _plot_bar_sisre_satellite_percentile(rpt, df, figure_dir, threshold=False, write_table=True, yrange=[0, 2])
        _plot_bar_sisre_satellite_percentile(rpt, df, figure_dir, threshold=True, write_table=False)
        # -TODO

        _add_to_report(rpt, figure_dir, df, df_month_perc, df_month_perc_rms, df_month_rms, file_vars)
        rpt.markdown_to_pdf()


def _add_to_report(
    rpt: "Report",
    figure_dir: PosixPath,
    df: pd.core.frame.DataFrame,
    df_month_perc: pd.core.frame.DataFrame,
    df_month_perc_rms: pd.core.frame.DataFrame,
    df_month_rms: pd.core.frame.DataFrame,
    file_vars: Dict[str, Any],
) -> None:
    """Add figures and tables to report

    Args:
        rpt:                Report object.
        figure_dir:          Figure directory.
        df:                 Given DAILY SISRE solutions are merged into one dataframe
        df_month_perc:      Dataframe with MONTHLY samples of 95th percentile SISRE
        df_month_perc_rms:  Dataframe with MONTHLY samples of 95th percentile SISRE, which are based on epochwise RMS
                            SISRE solutions
        df_month_rms:       Dataframe with MONTHLY samples of RMS SISRE
        file_vars:           File variables used for file and plot title naming.
    """

    rpt.add_text("\n\n##Monthly 95th percentile and RMS SISE for signal combinations (users)\n")

    rpt.add_text("95th percentile SISE results for signal combinations in meter: ")
    rpt.write_dataframe_to_markdown(df_month_perc, format="6.2f")
    rpt.add_figure(
        f"{figure_dir}/plot_bar_sisre_signal_combination_percentile.{FIGURE_FORMAT}\n",
        caption="Monthly 95th percentile of global average SISE for single- and dual-frequency users",
        clearpage=True,
    )
    rpt.add_figure(
        f"{figure_dir}/plot_bar_sisre_signal_combination_percentile_threshold.{FIGURE_FORMAT}\n",
        caption="Monthly 95th percentile of global average SISE for single- and dual-frequency users with red MPL threshold line",
        clearpage=True,
    )

    rpt.add_text(
        "95th percentile SISE results for signal combinations in meter based on epochwise RMS SISRE solutions: "
    )
    rpt.write_dataframe_to_markdown(df_month_perc_rms, format="6.2f")
    rpt.add_figure(
        f"{figure_dir}/plot_bar_sisre_signal_combination_percentile_rms.{FIGURE_FORMAT}\n",
        caption="Monthly 95th percentile of global average SISE for single- and dual-frequency users based on epochwise RMS SISRE solutions",
        clearpage=True,
    )
    rpt.add_figure(
        f"{figure_dir}/plot_bar_sisre_signal_combination_percentile_rms_threshold.{FIGURE_FORMAT}\n",
        caption="Monthly 95th percentile of global average SISE for single- and dual-frequency users based on epochwise RMS SISRE solutions with red MPL threshold line",
        clearpage=True,
    )

    rpt.add_text("RMS SISE results for signal combinations in meter: ")
    rpt.write_dataframe_to_markdown(df_month_rms, format="6.2f")
    rpt.add_figure(
        f"{figure_dir}/plot_bar_sisre_signal_combination_rms.{FIGURE_FORMAT}\n",
        caption="Monthly RMS of global average SISE for single- and dual-frequency users",
        clearpage=True,
    )
    rpt.add_text("\\newpage\n")


def _generate_dataframe(dsets: Dict[str, "Dataset"], field: str) -> Tuple[pd.core.frame.DataFrame]:
    """Generate dataframes based on given FIELD of datasets

    The dataframe "df" has following columns:

        time_gps:       Time in GPS time scale given as datetime objects
        satellite:      Satellite identifiers
        system:         GNSS identifier
        <solution_1>:   First FIELD solution (e.g. E1)
        <solution_2>:   Second FIELD solution (e.g. E1/E5b)
        <solution_3>:   Second FIELD solution (e.g. E1/E5a)

    Example for "df" dictionary:
     
                           time_gps satellite system        E1    E1/E5b    E1/E5a
        0       2019-01-01 00:00:00       E01      E  0.173793  0.123220  0.171849
        1       2019-01-01 00:00:00       E02      E  0.048395  0.127028  0.108108
        2       2019-01-01 00:00:00       E03      E  0.089328  0.121884  0.079576
        3       2019-01-01 00:00:00       E04      E  0.110866  0.088446  0.092292
        4       2019-01-01 00:00:00       E05      E  0.348935  0.305333  0.258733


    "df_month_perc" is a dataframe with month as indices and FIELD 95% percentile values for each signal combination
     as columns.

    Example for "df_month_perc" dictionary:

                        E1    E1/E5b    E1/E5a
        Jan-2019  0.335688  0.297593  0.326859
        Feb-2019  0.380575  0.330701  0.352535
        Mar-2019  0.353586  0.314817  0.344597

    The dataframe struture for "df_month_perc_rms" and "df_month_rms" is the same as for "df_month_perc".

    Args:
        dsets: Dictionary with FIELD solution name as keys (e.g. cnes_inav_e1, cnes_inav_e1e5b, cnes_fnav_e1e5a) and
               the belonging Dataset as value
        field: Field for which dataframe should be generated

    Returns:
        Tuple with following entries:

        | Element              | Description                                                                          |
        |----------------------|--------------------------------------------------------------------------------------|
        | df                   | Given DAILY FIELD solutions are merged into one dataframe                            |
        | df_month_perc        | Dataframe with MONTHLY samples of 95th percentile FIELD (based on Galileo SDD v1.0   |
        |                      | version)                                                                             |
        | df_month_perc_rms    | Dataframe with MONTHLY samples of 95th percentile FIELD, which are based on epochwise|
        |                      | RMS FIELD solutions (based on Galileo SDD v1.1 version)                              |
        | df_month_rms         | Dataframe with MONTHLY samples of RMS FIELD                                          |

    """
    user_types = config.tech[_SECTION].user_types.str
    user_types = user_types.split(",") if config.tech[_SECTION].user_types.str else user_types
    df = pd.DataFrame()
    column_names = list()
    df_month_perc_rms = None
    df_month_rms = None

    for idx, (name, dset) in enumerate(dsets.items()):
        
        if len(dset.meta["systems"]) > 1:
            log.fatal(
                f"The writer '{FILE_NAME}' can only be used, if each dataset is based only on one GNSS "
                f"(not '{', '.join(dset.meta['systems'])})'."
            )

        if dset.num_obs == 0:
            log.warn(f"Dataset '{name}' is empty.")
            continue

        # Get column names
        if user_types:
            column_name = user_types[idx]
        else:
            column_name = ""
            if "has" in dset.vars["id"] or "HAS" in dset.vars["id"]:
                column_name = "HAS" 
            elif "os" in dset.vars["id"] or "OS" in dset.vars["id"]:
                column_name = "OS"
            column_name = f"{column_name} {_get_signal_type(dset.meta)}" if column_name else _get_signal_type(dset.meta)
            
        column_names.append(column_name)
            
        df_tmp = dset.as_dataframe(fields=["satellite", "system", field, "time.gps"]) 
        df_tmp = df_tmp.rename(columns={field: column_name})

        if df.empty:
            df = df_tmp
            continue
        df = df.merge(df_tmp, on=["satellite", "system", "time_gps"], how="outer")

    if df.empty:
        log.fatal(f"All given datasets are empty [{', '.join(dsets.keys())}].")

    # Generate monthly samples of 95th percentile FIELD (after SDD v1.0 version)
    df_month_perc = df.drop(columns=["satellite", "system"]).set_index("time_gps").resample("M").apply(lambda x: np.nanpercentile(x, q=95))
    df_month_perc.index = df_month_perc.index.strftime("%b-%Y")

    # Generate monthly samples of RMS FIELD
    df_month_rms = df.drop(columns=["satellite", "system"]).set_index("time_gps").resample("M").apply(lambda x: np.sqrt(np.nanmean(np.square(x))))
    df_month_rms.index = df_month_rms.index.strftime("%b-%Y")

    # Generate monthly samples of 95th percentile FIELD based on epochwise FIELD RMS solutions(after SDD v1.1 version)
    #
    # NOTE: Following solutions assumes that each FIELD "column"-solution in dataframe 'df' is only given for one GNSS
    epochs = sorted(set(df["time_gps"]))
    df_tmp = pd.DataFrame(index=epochs, columns=column_names)

    # Loop over observation epochs
    for epoch in epochs:
        idx = df["time_gps"] == epoch
        row = dict()

        # Determine RMS for each signal type over all given FIELD satellite solutions in each epoch
        for name in column_names:
            row[name] = np.sqrt(np.nanmean(np.square(df[name][idx])))
        df_tmp.loc[epoch] = pd.Series(row)

    df_month_perc_rms = df_tmp.resample("M").apply(lambda x: np.nanpercentile(list(x), q=95))
    df_month_perc_rms.index = df_month_perc_rms.index.strftime("%b-%Y")
    df_month_perc_rms = df_month_perc_rms.transpose()

    return df, df_month_perc.transpose(), df_month_perc_rms, df_month_rms.transpose()


def _get_signal_type(meta: Dict[str, Any]) -> Tuple[str]:
    """Get signal type used for SISRE dataset

    Args:
        meta:   Dataset meta dictionary

    Returns:
        Signal type (e.g. Galileo E1, Galileo E1/E5a, GPS L1, ...)
    """
    system = meta["systems"][0]
    try:
        signal_type = meta["frequencies"][system].replace("_", "/")

    except KeyError:
        log.fatal(f"No frequencies are defined for GNSS {system!r} for option 'frequencies'.")

    return f"{enums.get_value('gnss_id_to_3digit_id', system)} {signal_type}"


#
# PLOT FUNCTIONS
#
def _plot_bar_sisre_signal_combination_percentile(
    df_month_perc: pd.core.frame.DataFrame, figure_dir: PosixPath, threshold: bool = False, suffix: str = ""
) -> None:
    """Generate bar plot with monthly SISRE 95% percentile for each GNSS signal combination

    Args:
       df_month_perc:  Dataframe with MONTHLY samples of 95th percentile SISRE
       figure_dir:      Figure directory.
       threshold:      Plot threshold.
       suffix:           File suffix.
    """
    df_month_perc.plot(kind="bar")
    # df_month_perc.plot(kind="bar", figsize=(5,3))

    if threshold:
        plt.axhline(2, color="r")
        plot_filename = f"plot_bar_sisre_signal_combination_percentile{suffix}_threshold.{FIGURE_FORMAT}"
    else:
        plot_filename = f"plot_bar_sisre_signal_combination_percentile{suffix}.{FIGURE_FORMAT}"
    plt.xlabel("Signal combination for single- and dual-frequency users")
    plt.xticks(rotation=0)
    plt.ylabel("SISE 95% [m]")
    # plt.legend(bbox_to_anchor=(1.04, 1), loc=2, borderaxespad=0., ncol=1) #legend on right side
    # plt.legend(bbox_to_anchor=(1, 1.15), loc=1, borderaxespad=0., ncol=3) #legend on top
    plt.legend(bbox_to_anchor=(0.8, -0.15), loc=1, borderaxespad=0.0, ncol=3)  # legend below
    # plt.legend(bbox_to_anchor=(0.8, -0.1), loc=1, borderaxespad=0., ncol=3) #legend below without xlabel
    plt.tight_layout()
    plt.savefig(figure_dir / plot_filename)
    plt.clf()  # clear the current figure


def _plot_bar_sisre_signal_combination_rms(df_month_rms: pd.core.frame.DataFrame, figure_dir: PosixPath):
    """Generate bar plot with monthly SISRE RMS for each GNSS signal combination

    Args:  
       df_month_rms:  Dataframe with MONTHLY samples of RMS SISRE
       figure_dir:     Figure directory.
    """
    df_month_rms.plot(kind="bar")

    plt.xlabel("Signal combination for single- and dual-frequency users")
    plt.xticks(rotation=0)
    plt.ylabel("SISE RMS [m]")
    # plt.legend(bbox_to_anchor=(1.04, 1), loc=2, borderaxespad=0., ncol=1) #legend on right side
    # plt.legend(bbox_to_anchor=(1, 1.15), loc=1, borderaxespad=0., ncol=3) #legend on top
    plt.legend(bbox_to_anchor=(0.8, -0.15), loc=1, borderaxespad=0.0, ncol=3)  # legend below
    # plt.legend(bbox_to_anchor=(0.8, -0.1), loc=1, borderaxespad=0., ncol=3) #legend below without xlabel
    plt.tight_layout()
    plt.savefig(figure_dir / f"plot_bar_sisre_signal_combination_rms.{FIGURE_FORMAT}")
    plt.clf()  # clear the current figure


def _plot_bar_sisre_satellite_percentile(
    rpt: "Report",
    df: pd.core.frame.DataFrame,
    figure_dir: PosixPath,
    threshold: bool = False,
    write_table: bool = False,
    yrange: Union[None, List[int]] = None,
) -> None:
    """Generate bar plot with monthly SISRE 95% percentile for each satellite

    Args:
       rpt:             Report object.
       df:              Dataframe with time, satellite, system and GNSS signal combinations as columns
       figure_dir:       Figure directory.
       threshold:       Plot threshold.
       write_table:     Write table.
       yrange:          Y-axis range of plot.
    """
    # Get user types by keeping order
    user_types = list(df.columns)
    for value in ["satellite", "system", "time_gps"]:
        user_types.remove(value)



    fig, axes = plt.subplots(len(user_types), 1, sharex=False, sharey=True)  # figsize=(6, 6));

    # Guarantee that 'axes' is iterable, which is needed by using 'zip' command.
    # Note: 'plt.subplots' does not return an iterable numpy array, if row and column number is one. This is the case,
    #       if 'user_types' has only one element.
    if not isinstance(axes, (np.ndarray)):
        axes = np.array([axes])

    if threshold:
        plot_filename = f"plot_bar_sisre_satellite_percentile_threshold.{FIGURE_FORMAT}"
    else:
        plot_filename = f"plot_bar_sisre_satellite_percentile.{FIGURE_FORMAT}"

    # plt.subplots_adjust(wspace=0.5, hspace=0.5);
    fig.set_figheight(10)  # inches

    # Generate subplots
    for ax, user in zip(axes, user_types):
        for sys in set(df.system):
            idx = df.system == sys
            if df[user][idx].isnull().all():
                log.debug(f"No valid data available for {enums.get_value('gnss_id_to_name', sys)} and column {user}.")
                continue
            df_user = df[idx].pivot(index="time_gps", columns="satellite", values=user)
            df_user_monthly_percentile = df_user.resample("M").apply(lambda x: np.nanpercentile(x, q=95))
            df_user_monthly_percentile.index = df_user_monthly_percentile.index.strftime("%b-%Y")
            df_user_monthly_percentile.transpose().plot(kind="bar", ax=ax, legend=False, title=user)
            ax.set_ylabel("SISE 95% [m]")
            if threshold:
                ax.axhline(7, color="r")
            if yrange is not None:
                ax.set_ylim(yrange)
    
            if write_table:
                rpt.add_text(f"95th percentile SISE results for signal combination **{user}** in meter for {enums.get_value('gnss_id_to_name', sys)}: ")
                rpt.write_dataframe_to_markdown(df_user_monthly_percentile.transpose(), format="6.2f")
                rpt.add_text("\\newpage\n")

    y_bbox = -0.5 if len(user_types) > 3 else -0.3
    plt.legend(bbox_to_anchor=(0.8, y_bbox), loc=1, borderaxespad=0.0, ncol=3)
    plt.tight_layout()
    plt.savefig(figure_dir / plot_filename)
    plt.clf()  # clear the current figure

    rpt.add_figure(
        f"{figure_dir}/{plot_filename}\n",
        caption="Monthly 95th percentile of global average SISE of each satellite for single- and dual-frequency users",
        clearpage=True,
    )

    rpt.add_text("\\newpage\n")
