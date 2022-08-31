"""Write report about a RINEX navigation file analysis run

Description:
------------

asdf


"""
# Standard library imports
from datetime import datetime
from typing import Tuple

# External library imports
import pandas as pd

# Midgard imports
from midgard.collections import enums
from midgard.dev import plugins
from midgard.plot.matplotext import MatPlotExt

# Where imports
from where.lib import config
from where.lib import log
from where.postprocessors.gnss_compare_tgd import gnss_compare_tgd
from where.writers._gnss_plot import GnssPlot
from where.writers._report import Report

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])

FIGURE_DPI = 200
FIGURE_FORMAT = "png"


@plugins.register
def rinex_nav_report(dset: "Dataset") -> None:
    """Write report about a RINEX navigation file analysis run

    Args:
        dset:        A dataset containing the data.
    """
    file_vars = {**dset.vars, **dset.analysis}
 
    # TODO: Better solution?
    if "station" not in file_vars:  # necessary if called for example by ./where/tools/concatenate.py
        file_vars["station"] = ""
        file_vars["STATION"] = ""

    # Generate figure directory to save figures generated for RINEX navigation file report
    figure_dir = config.files.path("output_rinex_nav_report_figure", file_vars=file_vars)
    figure_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate RINEX navigation file report
    path = config.files.path("output_rinex_nav_report", file_vars=file_vars)
    with config.files.open_path(path, create_dirs=True, mode="wt") as fid:
        rpt = Report(fid, rundate=dset.analysis["rundate"], path=path, description="RINEX navigation file analysis")
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
    # SIS status
    #
    if "sis_status" not in config.tech[_SECTION].skip_sections.list:
        rpt.add_text("\n# GNSS signal-in-space (SIS) status\n\n")
        plt = GnssPlot(dset, figure_dir, figure_format=FIGURE_FORMAT)

        # Plot GNSS SIS status (except if only Galileo system is available)
        if not (len(dset.unique("system")) == 1 and dset.unique("system")[0] == "E"):
            caption = "GNSS signal-in-space (SIS) status for each navigation message."
            if "E" in dset.unique("system"):  # Add extra comment for Galileo
                signal, nav_type = plt.get_first_galileo_signal()
                caption += (
                    f" Galileo SIS status is given for signal '{signal.upper()}' and "
                    f"navigation message type '{nav_type}'."
                )
            rpt.add_figure(
                    figure_path=plt.plot_gnss_signal_in_space_status_overview(), 
                    caption=caption, 
                    clearpage=True,
            )

        # Plots are only generated for Galileo
        for figure_path in plt.plot_gnss_signal_in_space_status():
            gnss = figure_path.stem.split("_")[5]
            if gnss == "galileo":
                caption=f"Galileo signal-in-space (SIS) status for signal {figure_path.stem.split('_')[-1].upper()}"
            else:
                caption=f"{gnss.upper()} signal-in-space (SIS) status"
                
            rpt.add_figure(
                    figure_path=figure_path, 
                    caption=caption, 
            )
            
        rpt.add_text("\n\\clearpage\n\n")

    #
    # Number of GNSS navigation messages
    #
    if "num_messages" not in config.tech[_SECTION].skip_sections.list:
        rpt.add_text("\n# Number of GNSS navigation messages\n\n")

        # Generate tables
        df = _table_navigation_message_overview(dset)
        rpt.write_dataframe_to_markdown(df, format="6.0f")

        # Plot GNSS navigation message overview (except if only Galileo system is available)
        if not (len(dset.unique("system")) == 1 and dset.unique("system")[0] == "E"):
            path = figure_dir / f"plot_bar_navigation_num_msg.{FIGURE_FORMAT}"
            plt = MatPlotExt()
            plt.plot_bar_dataframe_columns(
                df,
                column="num_msg",
                path=path,
                xlabel="Satellite",
                ylabel="Number of navigation messages",
                options={
                    "colormap": "viridis",
                    "figsize": (15, 10), 
                    "fontsize": 16,
                    "plot_to": "file",
                },
            )
            rpt.add_figure(path, "Number of navigation messages for each GNSS satellite.")

        # Plot Galileo navigation message overviews
        if "E" in dset.unique("system"):
            ylabel = {
                "fnav": "Number of F/NAV navigation messages",
                "inav": "Number of I/NAV navigation messages",
                "inav_e1": "Number of I/NAV navigation messages\ntransmitted by E1 signal",
                "inav_e5b": "Number of I/NAV navigation messages\ntransmitted by E5b signal",
                "inav_e1e5b": "Number of I/NAV navigation messages\ntransmitted by E1/E5b signal",
            }

            df_galileo = df[df["label"] == "Galileo"]  # Keep only Galileo observations
            df_galileo = df_galileo.drop(columns=["label"])  # Drop label column -> otherwise legend is plotted

            for column in ylabel.keys():
                if sum(df[column]) == 0:  # Skip plotting
                    continue
                path = figure_dir / f"plot_bar_navigation_{column}.{FIGURE_FORMAT}"
                plt = MatPlotExt()
                plt.plot_bar_dataframe_columns(
                    df_galileo, 
                    column=column, 
                    path=path, 
                    xlabel="Satellite", 
                    ylabel=ylabel[column],
                    options={
                        "figsize": (7, 4), 
                        "plot_to": "file",
                    },
                )
                rpt.add_figure(path, f"{ylabel[column]} for each satellite.")

        rpt.add_text("\n\\clearpage\n\n")
  
  
    #
    # TGD/BGD comparison
    #
    if "tgd_bgd_comparison" not in config.tech[_SECTION].skip_sections.list:
        rpt.add_text("\n# TGD/BGD comparison\n\n")

        bias_comp_def = {"bgd_e1_e5a_diff", "bgd_e1_e5b_diff", "tgd_diff", "tgd_b1_b2_diff", "tgd_b1_b3_diff"}
        
        # Add DCB comparison results to dataset if not existing
        if not set(dset.fields).intersection(bias_comp_def):
            gnss_compare_tgd(dset)

        if set(dset.fields).intersection(bias_comp_def):
            plt = GnssPlot(dset, figure_dir, figure_format=FIGURE_FORMAT)
            for figure_path in plt.plot_tgd_comparison():
                words = figure_path.stem.split("_")
                gnss = words[2] if words[1] == "field" else words[1]
                
                if "diff." in str(figure_path):
                    caption=f"TGD/BGD comparison against DCBs for {enums.gnss_id_to_name[gnss].value}"
                elif "diff_mean." in str(figure_path):
                    caption=f"TGD/BGD comparison against DCBs for {enums.gnss_id_to_name[gnss].value} (zero mean)"
                else:
                    caption=f"TGD/BGD for {enums.gnss_id_to_name[gnss].value}"
                rpt.add_figure(
                        figure_path=figure_path, 
                        caption=caption, 
                )
            
        else:
            log.warn(f"No TGD/BGD comparison plots are generated.")




#
# TABLE GENERATION FUNCTIONS
#
def _table_navigation_message_overview(dset: "Dataset"):
    """Generate Dataframe table with overview over number of navigation messages

    Args:
        dset:      A dataset containing the data.

    Returns:
        Dataframe with satellites as indices and following columns:

        | Name        | Description                                                                                  |
        |-------------|----------------------------------------------------------------------------------------------|
        | num_msg     | Number of navigation messages for each satellite                                             |
        | fnav        | Number of Galileo F/NAV navigation messages for each satellite                               |
        | inav        | Number of Galileo F/NAV navigation messages for each satellite                               |
        | inav_e1     | Number of Galileo I/NAV navigation messages transmitted by E1 signal for each satellite      |
        | inav_e5b    | Number of Galileo I/NAV navigation messages transmitted by E5b signal for each satellite     |
        | inav_e1e5b  | Number of Galileo merged I/NAV navigation messages transmitted by E1 and E5b signal for each |
        |             | satellite                                                                                    |

        Example:

            |    |num_msg| fnav| inav| inav_e1| inav_e5b| inav_e1e5b|
            |----|-------|-----|-----|--------|---------|-----------|
            | C01|     72|    0|    0|       0|        0|          0|
            | C02|     72|    0|    0|       0|        0|          0|
            | .. |    ...|  ...|  ...|     ...|      ...|        ...|
            | E01|    637|  159|  478|     160|      159|        159|
            | E02|    615|  155|  460|     155|      150|        155|
            | .. |    ...|  ...|  ...|     ...|      ...|        ...|

    """
    df = dset.as_dataframe()

    # Keep only data for run date
    day_start, day_end = _get_day_limits(dset)
    df = df[(df.time >= day_start) & (df.time <= day_end)]

    # Generate dataframe table with overview over number of navigation messages
    columns = ["label", "num_msg", "fnav", "inav", "inav_e1", "inav_e5b", "inav_e1e5b"]
    df_nav = pd.DataFrame(columns=columns)
    for satellite in sorted(dset.unique("satellite")):

        label = enums.gnss_id_to_name[satellite[0]].value  # GNSS identifier used as label
        num_msg = len(df.query(f"satellite == '{satellite}'"))
        fnav = len(df.query(f"satellite == '{satellite}' and nav_type == 'FNAV_E5a'"))
        inav_e1 = len(df.query(f"satellite == '{satellite}' and nav_type == 'INAV_E1'"))
        inav_e5b = len(df.query(f"satellite == '{satellite}' and nav_type == 'INAV_E5b'"))
        inav_e1e5b = len(df.query(f"satellite == '{satellite}' and nav_type == 'INAV_E1E5b'"))
        inav = inav_e1 + inav_e5b + inav_e1e5b

        row = [label, num_msg, fnav, inav, inav_e1, inav_e5b, inav_e1e5b]
        df_nav = df_nav.append(pd.DataFrame([row], columns=columns, index=[satellite]))

    return df_nav


#        
# AUXILIARY FUNCTIONS
# 
def _get_day_limits(dset: "Dataset") -> Tuple[datetime, datetime]:
    """Get start and end time for given run date

        Args:
            dset:      A dataset containing the data.

        Returns:
            Start and end date. 
        """
    day_start = min(dset.time.datetime)
    day_end = max(dset.time.datetime)

    return day_start, day_end
