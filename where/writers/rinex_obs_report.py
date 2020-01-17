"""Write report about a RINEX observation file analysis run

Description:
------------

asdf


"""
# Standard library imports
from typing import Any, Dict, List, Tuple, Union

# External library imports
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
def rinex_obs_report(dset: "Dataset") -> None:
    """Write report about a RINEX observation file analysis run

    Args:
        dset:        A dataset containing the data.
    """

    # TODO: Better solution?
    if "station" not in dset.vars:  # necessary if called for example by ./where/tools/concatenate.py
        dset.vars["station"] = ""
        dset.vars["STATION"] = ""

    # Generate figure directory to save figures generated for RINEX observation file report
    figure_dir = config.files.path("output_rinex_obs_report_figure", file_vars=dset.vars)
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    _plot_scatter_satellite_availability(dset, figure_dir)
    _table_observation_overview(dset)

    # Generate RINEX observation file report
    path = config.files.path("output_rinex_obs_report", file_vars=dset.vars)
    with config.files.open_path(path, create_dirs=True, mode="wt") as fid:
        rpt = Report(fid, rundate=dset.analysis["rundate"], path=path, description="RINEX observation file analysis")
        rpt.title_page()
        rpt.write_config()
        # TODO _add_figures(dset, rpt, figure_dir)
        _add_tables(dset, rpt, figure_dir)
        rpt.markdown_to_pdf()


def _add_figures(dset: "Dataset", rpt: "Report", figure_dir: "pathlib.PosixPath") -> None:
    """Add figures to report

    Args:
        dset:        A dataset containing the data.
        rpt:         Report object.
        figure_dir:  Figure directory.
    """
    rpt.add_text("\n# GNSS satellite availability\n\n")

    # Plot satellite availability
    if not (len(dset.unique("system")) == 1 and dset.unique("system")[0] == "E"):
        figure_name = f"plot_satellite_availability.{FIGURE_FORMAT}"
        caption = "GNSS satellite availability."
        rpt.add_figure(f"{figure_dir}/{figure_name}", caption, clearpage=True)


def _add_tables(dset: "Dataset", rpt: "Report", figure_dir: "pathlib.PosixPath") -> None:
    """Add tables and related plots to report

    Args:
        dset:        A dataset containing the data.
        rpt:         Report object.
        figure_dir:  Figure directory.
    """
    # Generate tables
    df_obstype, df_system = _table_observation_overview(dset)

    # Generate GNSS observation overview by system
    rpt.add_text("\n# GNSS observation overview by system\n\n")
    rpt.write_dataframe_to_markdown(df_system, format="6.0f")

    path = figure_dir / f"plot_gnss_number_of_satellites.{FIGURE_FORMAT}"
    plot_bar_dataframe_columns(
        df_system,
        column="num_sat",
        path=path,
        xlabel="GNSS name",
        ylabel="Number of satellites",
        label="sys",
        opt_args={"legend": False},
    )
    rpt.add_figure(path, "Number of satellites for each GNSS.")

    path = figure_dir / f"plot_gnss_observation_overview.{FIGURE_FORMAT}"
    plot_bar_dataframe_columns(
        df_system,
        column="num_obs",
        path=path,
        xlabel="GNSS name",
        ylabel="Number of observations",
        label="sys",
        opt_args={"legend": False},
    )
    rpt.add_figure(path, "Number of observations for each GNSS.", clearpage=True)

    # Generate GNSS observation overview by observation type
    rpt.add_text("\n# GNSS observation overview by observation type\n\n")
    rpt.write_dataframe_to_markdown(df_obstype, format="6.0f")

    path = figure_dir / f"plot_gnss_obstype_number_of_satellites.{FIGURE_FORMAT}"
    plot_bar_dataframe_columns(
        df_obstype,
        column="num_sat",
        path=path,
        xlabel="Observation type",
        ylabel="Number of satellites",
        label="sys",
        opt_args={"figsize": (15, 10), "fontsize": 16},
    )
    rpt.add_figure(path, "Number of satellites for each GNSS observation type.")

    path = figure_dir / f"plot_gnss_observation_type_overview.{FIGURE_FORMAT}"
    plot_bar_dataframe_columns(
        df_obstype,
        column="num_obs",
        path=path,
        xlabel="Observation type",
        ylabel="Number of observations",
        label="sys",
        opt_args={"figsize": (15, 10), "fontsize": 16},
    )
    rpt.add_figure(path, "Number of observations for each GNSS observation type.", clearpage=True)


#
# TABLE GENERATION FUNCTIONS
#
def _table_observation_overview(dset: "Dataset") -> Tuple[pd.core.frame.DataFrame, pd.core.frame.DataFrame]:
    """Generate Dataframe table with overview over number of observations

    Args:
        dset:      A dataset containing the data.

    Returns:
        Tuple with following dataframes:

        | Name        | Type       | Description                                                                     |
        |-------------|----------------------------------------------------------------------------------------------|
        | df_obstype  | np.ndarray | GNSS observation overview in relation to observation type. Indices are          |
        |             |            | combination of GNSS identifier and observation type, e.g. 'G C1C' or 'E L8X'    |
        | df_system   | np.ndarray | GNSS observation overview in relation to GNSS. Indices are GNSS identifiers,    |
        |             |            | e.g. E, G, ..., osv.                                                            |
        
        whereby following columns are defined:

        | Name        | Description                                                                                  |
        |-------------|----------------------------------------------------------------------------------------------|
        | sys         | GNSS identifier                                                                              |
        | num_sat     | Number of satellites                                                                         |
        | num_obs     | Number of observations                                                                       |

        Example for df_obstype:

            | index | sys | num_sat | num_obs |
            |-------|-----|---------|---------|
            | C C2X |  C  |      23 |    2469 |
            | C C7X |  C  |      11 |    1114 |
            | C D2X |  C  |      23 |    1200 |
            |   ... | ... |     ... |     ... |
            | E C1X |  E  |      21 |    2151 |
            | E C5X |  E  |      20 |    2066 |
            |   ... | ... |     ... |     ... |

        and for df_system:

            | index | sys | num_sat | num_obs |
            |-------|-----|---------|---------|
            |     C |  C  |      23 |   11949 |
            |     E |  E  |      21 |   26072 | 
            |     G |  G  |      31 |   33958 |
            |     R |  R  |      23 |   30391 |
      
    """
    df = dset.as_dataframe()
    columns = ["sys", "num_sat", "num_obs"]
    df_obstype = pd.DataFrame(columns=columns)
    df_system = pd.DataFrame(columns=columns)

    # Generate dataframe table with overview over number of observations related to observation types and GNSS
    for system in sorted(dset.unique("system")):
        gnss_name = enums.gnss_id_to_name[system].value
        idx = dset.filter(system=system)

        # Generate observation overview related to observation types
        for obstype in _sort_string_array(dset.meta["obstypes"][system]):
            index = obstype
            idx_num = dset.obs[obstype][idx] != 0  # TODO: RINEX Parser should be changed to 'nan' instead!!!!
            num_sat = len(set(dset.satellite[idx][idx_num]))

            row = [
                gnss_name,  # GNSS name
                num_sat,  # Number of satellites
                len(dset.obs[obstype][idx][idx_num]),  # Number of observations
            ]
            df_obstype = df_obstype.append(pd.DataFrame([row], columns=columns, index=[index]))

        # Generate observation overview related to GNSS
        idx_sys = df_obstype["sys"] == gnss_name
        row = [
            gnss_name,  # GNSS name
            len(set(dset.satellite[idx])),  # Number of satellites
            df_obstype["num_obs"][idx_sys].sum(),  # Number of observations
        ]
        df_system = df_system.append(pd.DataFrame([row], columns=columns, index=[gnss_name]))

    return df_obstype, df_system


#
# PLOT FUNCTIONS
#
def _plot_scatter_satellite_availability(dset: "Dataset", figure_dir: "pathlib.PosixPath") -> None:
    """Generate GNSS satellite observation availability overview based on RINEX observation file

    Args:
       dset:        A dataset containing the data.
       figure_dir:  Figure directory
    """

    # Generate x- and y-axis data per satellite
    x_arrays = []
    y_arrays = []
    labels = []
    for satellite in sorted(dset.unique("satellite"), reverse=True):
        idx = dset.filter(satellite=satellite)
        x_arrays.append(dset.time.gps.datetime[idx])
        y_arrays.append(dset.satellite[idx])
        labels.append(set(dset.system[idx]).pop())

    # Plot scatter plot
    plot(
        x_arrays=x_arrays,
        y_arrays=y_arrays,
        xlabel="",
        ylabel="",
        y_unit="",
        labels=labels,
        figure_path=figure_dir / f"plot_satellite_availability.{FIGURE_FORMAT}",
        opt_args={
            "colormap": "tab20",
            "figsize": (6.5, 6.5),
            "legend": True,
            "plot_to": "file",
            "plot_type": "scatter",
            "title": "Satellite availability",
        },
    )


#
# AUXILIARY FUNCTIONS
#
def _sort_string_array(array: List[str]) -> List[str]:
    """Sort string array based on last two characters

    Args:
        array: String array
        
    Returns:
        Sorted string array
    """
    array.sort(key=lambda x: x[-2:3])
    return array
