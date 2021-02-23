"""Write report about a RINEX navigation file analysis run

Description:
------------

asdf


"""
# Standard library imports
from datetime import datetime, timedelta
from typing import Any, Dict, List, Tuple, Union

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.collections import enums
from midgard.dev import plugins
from midgard.plot.matplotlib_extension import plot_bar_dataframe_columns, plot

# Where imports
from where.lib import config
from where.writers._report import Report


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

    # Generate plots
    _plot_gnss_signal_in_space_status(dset, figure_dir)
    _plot_galileo_signal_in_space_status(dset, figure_dir)

    # Generate RINEX navigation file report
    path = config.files.path("output_rinex_nav_report", file_vars=file_vars)
    with config.files.open_path(path, create_dirs=True, mode="wt") as fid:
        rpt = Report(fid, rundate=dset.analysis["rundate"], path=path, description="RINEX navigation file analysis")
        rpt.title_page()
        rpt.write_config()
        _add_figures(dset, rpt, figure_dir)
        _add_tables(dset, rpt, figure_dir)
        rpt.markdown_to_pdf()


def _add_figures(dset: "Dataset", rpt: "Report", figure_dir: "pathlib.PosixPath") -> None:
    """Add figures to report

    Args:
        dset:        A dataset containing the data.
        rpt:         Report object.
        figure_dir:  Figure directory.
    """
    rpt.add_text("\n# GNSS signal-in-space (SIS) status\n\n")

    # Plot GNSS SIS status (except if only Galileo system is available)
    if not (len(dset.unique("system")) == 1 and dset.unique("system")[0] == "E"):
        figure_name = f"plot_gnss_signal_in_space_status.{FIGURE_FORMAT}"
        caption = "GNSS signal-in-space (SIS) status for each navigation message."
        if "E" in dset.unique("system"):  # Add extra comment for Galileo
            signal, nav_type = _get_first_galileo_signal(dset)
            caption += (
                f" Galileo SIS status is given for signal '{signal.upper()}' and "
                f"navigation message type '{nav_type}'."
            )
        rpt.add_figure(f"{figure_dir}/{figure_name}", caption, clearpage=True)

    # Plots are only generated for Galileo
    if "E" in dset.unique("system"):
        for signal in sorted(_select_galileo_signal(dset).keys()):
            figure_name = f"plot_galileo_signal_in_space_status_{signal}.{FIGURE_FORMAT}"

            caption = f"Galileo signal-in-space (SIS) status for signal {signal.upper()}."
            rpt.add_figure(f"{figure_dir}/{figure_name}", caption)

        rpt.add_text("\n\\clearpage\n\n")


def _add_tables(dset: "Dataset", rpt: "Report", figure_dir: "pathlib.PosixPath") -> None:
    """Add tables and related plots to report

    Args:
        dset:        A dataset containing the data.
        rpt:         Report object.
        figure_dir:  Figure directory.
    """

    rpt.add_text("\n# Number of GNSS navigation messages\n\n")

    # Generate tables
    df = _table_navigation_message_overview(dset)
    rpt.write_dataframe_to_markdown(df, format="6.0f")

    # Plot GNSS navigation message overview (except if only Galileo system is available)
    if not (len(dset.unique("system")) == 1 and dset.unique("system")[0] == "E"):
        path = figure_dir / f"plot_bar_navigation_num_msg.{FIGURE_FORMAT}"
        plot_bar_dataframe_columns(
            df,
            column="num_msg",
            path=path,
            xlabel="Satellite",
            ylabel="Number of navigation messages",
            opt_args={"colormap": "viridis", "figsize": (15, 10), "fontsize": 16},
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
            plot_bar_dataframe_columns(df_galileo, column=column, path=path, xlabel="Satellite", ylabel=ylabel[column])
            rpt.add_figure(path, f"{ylabel[column]} for each satellite.")


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
# PLOT FUNCTIONS
#
def _plot_galileo_signal_in_space_status(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Generate Galileo Signal-in-Space (SIS) status plot based on Galileo signal health status (SHS), SIS Accuracy 
    (SISA) and data validity status (DVS) given in RINEX navigation file.

    The SIS status can be:

     | CODE | SIS STATUS      | PLOTTED COLOR | DESCRIPTION                     |
     |------|-----------------|---------------|---------------------------------|
     |   0  | healthy         |         green | SIS status used by all GNSS     |
     |   1  | marginal (SISA) |        yellow | SIS status only used by Galileo |
     |   2  | marignal (DVS)  |        orange | SIS status only used by Galileo |
     |   3  | unhealthy       |           red | SIS status used by all GNSS     |

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory.
    """
    colors = ["green", "yellow", "orange", "red"]
    labels = ["healthy", "marginal (sisa)", "marginal (dvs)", "unhealthy"]
    status_def = [0, 1, 2, 3]
    signals = _select_galileo_signal(dset)

    # Generate plot for each given Galileo signal (e.g. E1, E5a, E5b)
    for signal, nav_type in sorted(signals.items()):

        x_arrays = []
        y_arrays = []
        for status in status_def:
            time, satellite = _get_gnss_signal_in_space_status_data(dset, status, signal, only_galileo=True)

            x_arrays.append(time)
            y_arrays.append(satellite)

        # Limit x-axis range to rundate
        day_start, day_end = _get_day_limits(dset)

        # Generate plot
        plot(
            x_arrays=x_arrays,
            y_arrays=y_arrays,
            xlabel="Time [GPS]",
            ylabel="Satellite",
            y_unit="",
            labels=labels,
            colors=colors,
            figure_path=figure_dir / f"plot_galileo_signal_in_space_status_{signal}.{FIGURE_FORMAT}",
            opt_args={
                "figsize": (7, 5),
                "marker": "s",
                "marksersize": 10,
                "legend_ncol": 4,
                "legend_location": "bottom",
                "plot_to": "file",
                "plot_type": "scatter",
                "title": f"Galileo signal-in-space status for signal {signal.upper()} ({nav_type})",
                "xlim": [day_start, day_end],
            },
        )


def _plot_gnss_signal_in_space_status(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Generate GNSS Signal-in-Space (SIS) status plot based on SIS status given in RINEX navigation file

    The SIS status can be:

     | CODE | SIS STATUS      | PLOTTED COLOR | DESCRIPTION                     |
     |------|-----------------|---------------|---------------------------------|
     |   0  | healthy         |         green | SIS status used by all GNSS     |
     |   1  | marginal (SISA) |        yellow | SIS status only used by Galileo |
     |   2  | marignal (DVS)  |        orange | SIS status only used by Galileo |
     |   3  | unhealthy       |           red | SIS status used by all GNSS     |

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """
    colors = ["green", "yellow", "orange", "red"]
    labels = ["healthy", "marginal (sisa)", "marginal (dvs)", "unhealthy"]
    status_def = [0, 1, 2, 3]
    signal = None

    # Select only one Galileo signal
    # Note: Navigation message for Galileo can include I/NAV and F/NAV messages for different signals (E1, E5a, E5b).
    #       For plotting we choose only one of them.
    if "E" in dset.unique("system"):
        signal, _ = _get_first_galileo_signal(dset)

    # Generate time and satellite data for given SIS status
    x_arrays = []
    y_arrays = []
    for status in status_def:
        time, satellite = _get_gnss_signal_in_space_status_data(dset, status, signal)
        x_arrays.append(time)
        y_arrays.append(satellite)

    # Limit x-axis range to rundate
    day_start, day_end = _get_day_limits(dset)

    # Generate plot
    plot(
        x_arrays=x_arrays,
        y_arrays=y_arrays,
        xlabel="Time [GPS]",
        ylabel="Satellite",
        y_unit="",
        labels=labels,
        colors=colors,
        figure_path=figure_dir / f"plot_gnss_signal_in_space_status.{FIGURE_FORMAT}",
        opt_args={
            "figsize": (7, 11),
            "marker": "s",
            "marksersize": 10,
            "legend_ncol": 4,
            "legend_location": "bottom",
            "plot_to": "file",
            "plot_type": "scatter",
            "tick_labelsize": ("y", 7),  # specify labelsize 7 for y-axis
            "title": f"GNSS signal-in-space status",
            "xlim": [day_start, day_end],
        },
    )
    # TODO: Legend has to be improved. Old configuration:
    # figsize = (7, 10)
    # loc="lower center",
    # bbox_to_anchor=(0.5, -0.01),
    # frameon=True,


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
    if config.where.rinex_nav_report.only_for_rundate.bool:
        day_start = datetime(
            dset.analysis["rundate"].year, dset.analysis["rundate"].month, dset.analysis["rundate"].day
        )
        day_end = day_start + timedelta(seconds=86399)
    else:
        day_start = min(dset.time.datetime)
        day_end = max(dset.time.datetime)

    return day_start, day_end


def _get_first_galileo_signal(dset):
    """Get first Galileo signal given in the navigation message based on ordered list

    Args:
       dset:        A dataset containing the data.

    Returns:
       Tuple with chosen Galileo signal and navigation message type
    """
    signals = _select_galileo_signal(dset)
    signal = next(iter(sorted(signals)))
    nav_type = signals[signal]

    return signal, nav_type


def _get_gnss_signal_in_space_status_data(
    dset: "Dataset", status: int, signal: str, only_galileo: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """Get GNSS signal in space status (SIS) time and satellite data based on a given SIS status

    The SIS status can be:

     | CODE | SIS STATUS      | 
     |------|-----------------|
     |   0  | healthy         |
     |   1  | marginal (SISA) | 
     |   2  | marignal (DVS)  |
     |   3  | unhealthy       |

    The SIS status code is used for representing the SIS status in the dataset fields "sis_status_<signal>" (with 
    <signal>: e1, e5a or e5b).

    Args:
       dset:          A dataset containing the data.
       status:        Signal in space status.   
       signal:        Galileo signal for which SIS status should be determined. The signal can be e1, e5a or e5b.
       only_galileo:  Generate SIS status data only for Galileo

    Returns:
        Tuple with time and satellite data for a given SIS status and ordered by satellite
    """
    # TODO: order of satellites is not correct. Maybe to save the data in a dataframe could help?

    # Generate x- and y-axis data
    time = []
    satellite = []

    for sat in sorted(dset.unique("satellite"), reverse=True):
        idx = dset.filter(satellite=sat)

        if sat.startswith("E"):
            idx_status = getattr(dset, "sis_status_" + signal)[idx] == status
        else:
            if not only_galileo:
                if status == 0:  # healthy status
                    idx_status = dset.sv_health[idx] == 0
                elif status == 3:  # unhealthy status
                    idx_status = dset.sv_health[idx] > 0
                else:
                    continue
            else:
                continue

        time.extend(dset.time.gps.datetime[idx][idx_status])
        satellite.extend(dset.satellite[idx][idx_status])

    return time, satellite


def _select_galileo_signal(dset: "Dataset") -> Tuple[List[str], str]:
    """Select Galileo signal depending on given data in Dataset
    
    Args:
        dset: A dataset containing the data.

    Returns:
        Selected Galileo signal
    """
    signals = dict()

    for nav_type in dset.unique("nav_type"):
        if nav_type == "FNAV_E5a":
            signals.update({"e5a": "F/NAV"})
            continue
        elif nav_type == "INAV_E1":
            signals.update({"e1": "I/NAV"})
            continue
        elif nav_type == "INAV_E5b":
            signals.update({"e5b": "I/NAV"})
            continue
        elif nav_type == "INAV_E1E5b":
            signals.update({"e1": "I/NAV", "e5b": "I/NAV"})
            continue

    return signals
