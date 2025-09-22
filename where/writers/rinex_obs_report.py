"""Write report about a RINEX observation file analysis run

Description:
------------

asdf


"""
# Standard library imports
import pathlib
from typing import List, Tuple

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.collections import enums
from midgard.dev import plugins
from midgard.plot.matplotext import MatPlotExt

# Where imports
from where.data import dataset3 as dataset
from where.lib import config, log
from where.writers._gnss_plot import GnssPlot
from where.writers._report import Report

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


FIGURE_DPI = 200
FIGURE_FORMAT = "png"


@plugins.register
def rinex_obs_report(dset: "Dataset") -> None:
    """Write report about a RINEX observation file analysis run

    Args:
        dset:        A dataset containing the data.
    """
    file_vars = {**dset.vars,**dset.analysis}

    # TODO: Better solution?
    if "station" not in file_vars:  # necessary if called for example by ./where/tools/concatenate.py
        file_vars["station"] = ""
        file_vars["STATION"] = ""

    # Generate figure directory to save figures generated for RINEX observation file report
    figure_dir = config.files.path("output_rinex_obs_report_figure", file_vars=file_vars)
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Get 'read' Dataset
    dset_read = dataset.Dataset.read(**dict(file_vars, stage="read"))

    # Generate plots based on 'read' stage data
    df_system, df_obstype = _generate_dataframes(dset_read)
    _plot_observation_system(df_system, figure_dir)
    _plot_observation_type(df_obstype, figure_dir)
    
    # Generate RINEX observation file report
    path = config.files.path("output_rinex_obs_report", file_vars=file_vars)
    log.info(f"Write file {path}.")
    with config.files.open_path(path, create_dirs=True, mode="wt") as fid:
        rpt = Report(fid, rundate=dset.analysis["rundate"], path=path, description="RINEX observation file analysis")
        rpt.title_page()
        rpt.write_config()
        _add_to_report(rpt, figure_dir, dset_read, dset, df_system, df_obstype)
        rpt.markdown_to_pdf()


def _add_to_report(
        rpt: "Report", 
        figure_dir: pathlib.PosixPath,
        dset_read: "Dataset",
        dset_edit: "Dataset",
        df_system: pd.core.frame.DataFrame,
        df_obstype: pd.core.frame.DataFrame,
) -> None:
    """Add figures and tables to report

    Args:
        rpt:         Report object.
        figure_dir:  Figure directory.
        dset_read:   A dataset containing the data from the 'read' stage.
        dset_edit:   A dataset containing the data from the 'edit' stage.
        df_system:   Dataframe with GNSS observation overview in relation to GNSS. Indices are GNSS identifiers, e.g. 
                     E, G, ..., osv.. 
        df_obstype:  Dataframe with GNSS observation overview in relation to observation type. Indices are combination
                     of GNSS identifier and observation type, e.g. 'G C1C' or 'E L8X'.
    """
    plt_read = GnssPlot(dset_read, figure_dir)
    plt_edit = GnssPlot(dset_edit, figure_dir)

    #
    # Observation overview
    #            
    if "observation_overview" in config.tech[_SECTION].add_to_report.list:

        rpt.add_text("\n# Observation overview\n\n")        
                
        # Generate GNSS observation overview by system
        rpt.add_text("\n## GNSS observation overview by system\n\n")
        rpt.write_dataframe_to_markdown(df_system, format="6.0f")
    
        rpt.add_figure(
                figure_path=figure_dir / f"plot_gnss_number_of_satellites.{FIGURE_FORMAT}", 
                caption="Number of satellites for each GNSS.",
        )
        
        rpt.add_figure(
                figure_path=figure_dir / f"plot_gnss_observation_overview.{FIGURE_FORMAT}", 
                caption="Number of observations for each GNSS.", 
                clearpage=True,
        )
    
        # Generate GNSS observation overview by observation type
        rpt.add_text("\n## GNSS observation overview by observation type\n\n")
        rpt.write_dataframe_to_markdown(df_obstype, format="6.0f")
    
    
        rpt.add_figure(
                figure_path=figure_dir / f"plot_gnss_obstype_number_of_satellites.{FIGURE_FORMAT}",
                caption="Number of satellites for each GNSS observation type.",
        )
    
        rpt.add_figure(
                figure_path=figure_dir / f"plot_gnss_obstype_overview.{FIGURE_FORMAT}", 
                caption="Number of observations for each GNSS observation type.", 
                clearpage=True,
        )
     
    if "obstype_availability" in config.tech[_SECTION].add_to_report.list:
        for figure_path in plt_read.plot_obstype_availability():
            rpt.add_figure(
                    figure_path=figure_path, 
                    caption=f"Observation type availability for {enums.gnss_id_to_name[figure_path.stem[-1]]}", 
                    clearpage=True,
            )
        
        


    #
    # Satellite plot
    #        
    if "satellite_plot" in config.tech[_SECTION].add_to_report.list:        
        rpt.add_text("\n# Satellite plots\n\n")
    
        # Generate satellite availability overview
        rpt.add_text("\n## Satellite availability\n\n")
        rpt.add_figure(
                figure_path=plt_read.plot_satellite_availability(), 
                caption="Satellite availability.",
                clearpage=False,
        )
        
        rpt.add_figure(
                figure_path=plt_read.plot_number_of_satellites(), 
                caption="Number of satellites for each observation epoch per system.",
                clearpage=False,
        )


        for figure_path in plt_edit.plot_skyplot():
            rpt.add_figure(
                    figure_path=figure_path, 
                    caption=f"Skyplot for {enums.gnss_id_to_name[figure_path.stem[-1]]}", 
                    clearpage=False,
            )
        
        for figure_path in plt_edit.plot_satellite_elevation():
            rpt.add_figure(
                    figure_path=figure_path,
                    caption=f"Satellite elevation for {enums.gnss_id_to_name[figure_path.stem[-1]]}", 
                    clearpage=True,
            )
    
    #
    # Linear combination 
    #
    if "gnss_linear_combination" in config.tech.postprocessors.list:
        if config.tech.gnss_linear_combination.linear_combination.list:
            rpt.add_text("\n# Linear combination\n\n") 
            for fig in plt_edit.plot_linear_combinations():
                
                lin_comb = fig.name.replace("_", "-").capitalize() 
                system = enums.gnss_id_to_name[fig.system]
                rpt.add_figure(
                        figure_path=fig.file_path,
                        caption=f"{lin_comb} linear combination for {system} observations ({', '.join(fig.obstype)})", 
                        clearpage=True,
                )
            
            
    #
    # Epoch by epoch difference
    #
    if "gnss_epoch_by_epoch_difference" in config.tech.postprocessors.list:
        rpt.add_text("\n# Epoch by epoch difference\n\n")          
        for figure_path in plt_edit.plot_epoch_by_epoch_difference():
            system, obstype = figure_path.stem.split("_")[5:7] 
            rpt.add_figure(
                figure_path=figure_path, 
                caption=f"Epoch by epoch difference for {enums.gnss_id_to_name[system].value} observation {obstype}", 
                clearpage=True,
            )

    #
    # Fields to plot
    #
    fields_to_plot = config.tech[_SECTION].fields_to_plot.list
    if fields_to_plot:

        for field in fields_to_plot:
            
            if "." in field:
                collection, field = field.split(".")
            else:
                collection = None
                    
            for figure_path in plt_edit.plot_field(field, collection):
                system, field = figure_path.stem.split("_")[2:4] 
                rpt.add_figure(
                    figure_path=figure_path, 
                    caption=f"Field {field} for {enums.gnss_id_to_name[system].value} observation", 
                    clearpage=True,
                )
        
        

#
# GENERATE DATAFRAMES
#
def _generate_dataframes(dset: "Dataset") -> Tuple[pd.core.frame.DataFrame, pd.core.frame.DataFrame]:
    """Generate Dataframe table with overview over number of observations

    Args:
        dset:      A dataset containing the data.

    Returns:
        Tuple with following dataframes:

        | Name        | Type       | Description                                                                     |
        |-------------|----------------------------------------------------------------------------------------------|
        | df_system   | np.ndarray | GNSS observation overview in relation to GNSS. Indices are GNSS identifiers,    |
        |             |            | e.g. E, G, ..., osv.                                                            |
        | df_obstype  | np.ndarray | GNSS observation overview in relation to observation type. Indices are          |
        |             |            | combination of GNSS identifier and observation type, e.g. 'G C1C' or 'E L8X'    |
        
        whereby following columns are defined:

        | Name        | Description                                                                                  |
        |-------------|----------------------------------------------------------------------------------------------|
        | sys         | GNSS identifier                                                                              |
        | num_sat     | Number of satellites                                                                         |
        | num_obs     | Number of observations                                                                       |

        Example for df_system:

            | index | sys | num_sat | num_obs |
            |-------|-----|---------|---------|
            |     C |  C  |      23 |   11949 |
            |     E |  E  |      21 |   26072 | 
            |     G |  G  |      31 |   33958 |
            |     R |  R  |      23 |   30391 |
            
        and for df_obstype:

            | index | sys | num_sat | num_obs |
            |-------|-----|---------|---------|
            | C C2X |  C  |      23 |    2469 |
            | C C7X |  C  |      11 |    1114 |
            | C D2X |  C  |      23 |    1200 |
            |   ... | ... |     ... |     ... |
            | E C1X |  E  |      21 |    2151 |
            | E C5X |  E  |      20 |    2066 |
            |   ... | ... |     ... |     ... |
      
    """
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
            idx_num = np.logical_not(np.isnan(dset.obs[obstype][idx])) 
            num_sat = len(set(dset.satellite[idx][idx_num]))

            row = [
                gnss_name,  # GNSS name
                num_sat,  # Number of satellites
                len(dset.obs[obstype][idx][idx_num]),  # Number of observations
            ]
            df_obstype = pd.concat([df_obstype, pd.DataFrame([row], columns=columns, index=[index])])

        # Generate observation overview related to GNSS
        idx_sys = df_obstype["sys"] == gnss_name
        row = [
            gnss_name,  # GNSS name
            len(set(dset.satellite[idx])),  # Number of satellites
            df_obstype["num_obs"][idx_sys].sum(),  # Number of observations
        ]
        df_system = pd.concat([df_system, pd.DataFrame([row], columns=columns, index=[gnss_name])])

    return df_system, df_obstype


#
# PLOT FUNCTIONS
#
def _plot_observation_system(
        df_system: pd.core.frame.DataFrame,
        figure_dir: pathlib.PosixPath,
) -> None:
    """Plot observation system plots

    Args:
       df_system:   Dataframe with GNSS observation overview in relation to GNSS. Indices are GNSS identifiers, e.g. 
                    E, G, ..., osv..                               
       figure_dir:  Figure directory
    """ 
     
    plt = MatPlotExt()
    plt.plot_bar_dataframe_columns(
        df_system,
        column="num_sat",
        path=figure_dir / f"plot_gnss_number_of_satellites.{FIGURE_FORMAT}",
        xlabel="GNSS name",
        ylabel="# satellites",
        label="sys",
        options={
              "figsize": (8,4),
              "legend": False,
              "plot_to": "file",
        },
    )
    
    plt = MatPlotExt()
    plt.plot_bar_dataframe_columns(
        df_system,
        column="num_obs",
        path=figure_dir / f"plot_gnss_observation_overview.{FIGURE_FORMAT}",
        xlabel="GNSS name",
        ylabel="# observations",
        label="sys",
        options={
              "figsize": (8,4),
              "legend": False,
              "plot_to": "file",
        },
    )
    
    
def _plot_observation_type(
        df_obstype: pd.core.frame.DataFrame,
        figure_dir: pathlib.PosixPath,
) -> None:
    """Plot observation type plots

    Args:
       df_obstype:  Dataframe with GNSS observation overview in relation to observation type. Indices are combination
                    of GNSS identifier and observation type, e.g. 'G C1C' or 'E L8X'.
       figure_dir:  Figure directory
    """   
    num_sys = len(set(df_obstype.sys))


    plt = MatPlotExt()
    plt.plot_bar_dataframe_columns(
        df_obstype,
        column="num_sat",
        path=figure_dir / f"plot_gnss_obstype_number_of_satellites.{FIGURE_FORMAT}",
        xlabel="Observation type",
        ylabel="# satellites",
        label="sys",
        options={
              "figsize": (8, 4), 
              "fontsize": 17,
              "legend": True,
              "legend_location": "bottom",
              "legend_ncol": num_sys,
              "plot_to": "file",
        },
    )
    
    plt = MatPlotExt()
    plt.plot_bar_dataframe_columns(
        df_obstype,
        column="num_obs",
        path=figure_dir / f"plot_gnss_obstype_overview.{FIGURE_FORMAT}",
        xlabel="Observation type",
        ylabel="# observations",
        label="sys",
        options={
              "figsize": (8, 4), 
              "fontsize": 17,
              "legend": True,
              "legend_location": "bottom",
              "legend_ncol": num_sys,
              "plot_to": "file",
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
