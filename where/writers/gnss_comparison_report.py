"""Compare different GNSS Where datasets and write a report

Description:
------------

A dictionary with datasets is used as input for this writer. The keys of the dictionary are station names. 

Example:
--------

    from where import data
    from where import writers

    # Read a dataset
    dset = data.Dataset(rundate=rundate, tech=tech, stage=stage, dataset_name=name, dataset_id=dataset_id)

    # Write dataset
    writers.write_one('gnss_comparison_report', dset=dset, do_report=False)

"""
# Standard library imports
import copy
from typing import Any, Dict
from pathlib import PosixPath

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.data import position
from midgard.dev import plugins
from midgard.plot.matplotext import MatPlotExt

# Where imports
from where.lib import config
from where.lib import log
from where.writers._report import Report

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])

FIGURE_FORMAT = "png"
FILE_NAME = __name__.split(".")[-1]


@plugins.register
def gnss_comparison_report(dset: Dict[str, "Dataset"]) -> None:
    """Compare GNSS datasets

    Args:
        dset:  Dictionary with station name as keys and the belonging Dataset as value
    """
    dset_first = dset[list(dset.keys())[0]]
    file_vars = {**dset_first.vars, **dset_first.analysis}
    file_vars["solution"] = config.tech[_SECTION].solution.str.lower()

    # Generate figure directory to save figures generated for GNSS report
    figure_dir = config.files.path("output_gnss_comparison_report_figure", file_vars=file_vars)
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    _, dfs_day, dfs_month = _generate_dataframes(dset)
    _plot_position_error(dfs_day, dfs_month, figure_dir, file_vars)

    # Generate GNSS comparison report
    path = config.files.path("output_gnss_comparison_report", file_vars=file_vars)
    with config.files.open_path(path, create_dirs=True, mode="wt") as fid:
        rpt = Report(fid, rundate=dset_first.analysis["rundate"], path=path, description="Comparison of GNSS analyses")
        rpt.title_page()
        _add_to_report(rpt, figure_dir, dfs_day, dfs_month, file_vars)
        rpt.markdown_to_pdf()


#
# AUXILIARY FUNCTIONS
#
def _apply(df: pd.core.frame.DataFrame, sample: str, func: str) -> pd.core.frame.DataFrame:
    """Resample dataframe and apply given function 

    Args:
        df:      Dataframe.
        sample:  Sample definition ("D": day, "M": month)
        func:    Function to be applied ("mean", "percentile_68", "percentile_90", "percentile_95", "rms", "std")

    Returns:
        Resampled dataframe by applying given function
    """
    df_sampled = df.set_index("time_gps").resample(sample)

    if func == "mean":
        df_sampled = df_sampled.mean()
    elif func == "percentile_68":
        df_sampled = df_sampled.apply(lambda x: np.nanpercentile(x, q=68))
    elif func == "percentile_90":
        df_sampled = df_sampled.apply(lambda x: np.nanpercentile(x, q=90))
    elif func == "percentile_95":
        df_sampled = df_sampled.apply(lambda x: np.nanpercentile(x, q=95))
    elif func == "rms":
        df_sampled = df_sampled.apply(lambda x: np.sqrt(np.nanmean(np.square(x))))
    elif func == "std":
        df_sampled = df_sampled.std()
    else:
        log.fatal(f"Function '{func}' is not defined.")

    return df_sampled  


def _add_to_report(
    rpt: "Report",
    figure_dir: PosixPath,
    dfs_day: Dict[str, pd.core.frame.DataFrame],
    dfs_month: Dict[str, pd.core.frame.DataFrame],
    file_vars: Dict[str, Any],
) -> None:
    """Add figures and tables to report

    Args:
        rpt:         Report object.
        figure_dir:  Figure directory.
        dfs_day:     Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a dictionary as
                     values. The dictionary has fields as keys (e.g. hpe, vpe) and the belonging dataframe as value with
                     DAILY samples of e.g. 95th percentile and stations as columns.
        dfs_month    Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a dictionary as
                     values. The dictionary has fields as keys (e.g. hpe, vpe) and the belonging dataframe as value with
                     MONTHLY samples of e.g. 95th percentile and stations as columns.
        file_vars:   File variables used for file and plot title naming.
    """

    text_def = {
            "mean": "average",
            "percentile_68": "68th percentile",
            "percentile_90": "90th percentile",
            "percentile_95": "95th percentile",
            "std": "standard deviation",
            "rms": "RMS",
    }

    field_def = {
            "east": "East component of topocentric coordinates",
            "north": "North component of topocentric coordinates",
            "up":  "Up component of topocentric coordinates",
            "hpe": "horizontal position error (HPE)",
            "vpe": "vertical position error (VPE)",
            "pos_3d": "3D position error (3D)",
            "pdop": "position dilution of precision (PDOP)",
            "hdop": "horizontal dilution of precision (HDOP)",
            "vdop": "vertical dilution of precision (VDOP)",
    }

    for type_ in dfs_day.keys():

        for sample in config.tech[_SECTION].samples.list:

            sample = sample.capitalize()
            rpt.add_text(f"\n# {sample} {text_def[type_]} for given solutions\n\n")
      
            if sample == "Daily":
                for field in config.tech[_SECTION].fields.list:
                    dfs_day[type_][field].index = dfs_day[type_][field].index.strftime("%d-%m-%Y")
                    rpt.add_text(f"Daily {text_def[type_]} {field.upper()} results in meter:")
                    rpt.write_dataframe_to_markdown(dfs_day[type_][field], format="6.2f", statistic=True)


            elif sample == "Monthly":
                for field in config.tech[_SECTION].fields.list:
                    rpt.add_text(f"Monthly {text_def[type_]} {field.upper()} results in meter:")
                    rpt.write_dataframe_to_markdown(dfs_month[type_][field], format="6.2f")


            # Add plots
            for field in config.tech[_SECTION].fields.list:
                rpt.add_figure(
                    f"{figure_dir}/plot_{type_}_{field}_{sample.lower()}_{file_vars['date']}_{file_vars['solution'].lower()}.{FIGURE_FORMAT}",
                    caption=f"{text_def[type_].capitalize()} for {field_def[field]}.",
                    clearpage=True,
                )


def _generate_dataframes(dset: Dict[str, "Dataset"]) -> Dict[str, pd.core.frame.DataFrame]:
    """Generate dataframe based on station datasets

    The dataframe for each station in dictionary "dfs" has following columns:

        east:   East-coordinate in topocentric system
        north:  North-coordinate in topocentric system
        up:     Up-coordinate in topocentric system
        hpe:    horizontal position error
        vpe:    vertical position error
        pos_3d: 3D position error
        pdop:   position dilution of precision
        hdop:   horizontal dilution of precision
        vdop:   vertical dilution of precision

    Example for "dfs" dictionary:
     
             'hons':                   time_gps       hpe       vpe      east     north        up
                        0      2019-03-01 00:00:00  0.301738  0.057244  0.113758  0.279472  0.057244
                        1      2019-03-01 00:00:00  0.301738  0.057244  0.113758  0.279472  0.057244

             'krss':                   time_gps       hpe       vpe      east     north        up
                        0      2019-03-01 00:00:00  0.710014  0.186791 -0.235267  0.669903  0.186791
                        1      2019-03-01 00:00:00  0.710014  0.186791 -0.235267  0.669903  0.186791

    Example for "dfs_day" dictionary for "mean" key:
        'mean':{

             'hpe':                 nabf      vegs      hons      krss
                        time_gps                                          
                        2019-03-01  1.368875  0.935687  1.136763  0.828754
                        2019-03-02  0.924839  0.728280  0.911677  0.854832


             'vpe':                 nabf      vegs      hons      krss
                        time_gps                                          
                        2019-03-01  1.715893  1.147265  1.600330  0.976541
                        2019-03-02  1.533437  1.307373  1.476295  1.136991
        }

    Example for "dfs_month" dictionary for "mean" key:
        'mean':{
            'hpe':                nabf      vegs      hons      krss
                        Mar-2019  1.186240  0.861718  1.095827  1.021354
                        Apr-2019  0.891947  0.850343  0.977908  0.971099

            'vpe':                nabf      vegs      hons      krss
                        Mar-2019  1.854684  1.291406  1.450466  1.225467
                        Apr-2019  1.964404  1.706507  1.687994  1.500742
        }


    Args:
        dset: Dictionary with station name as keys and the belonging Dataset as value

    Returns:
        Tuple with following entries:

        | Element              | Description                                                                          |
        |----------------------|--------------------------------------------------------------------------------------|
        | dfs                  | Dictionary with station name as keys and the belonging dataframe as value with       |
        |                      | following dataframe columns: east, north, up, hpe, vpe, pos_3d                       |
        | dfs_day              | Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a  |
        |                      | dictionary as values. The dictionary has fields as keys (e.g. hpe, vpe) and the      |
        |                      | belonging dataframe as value with DAILY samples of e.g. 95th percentile_xx and       | 
        |                      | stations as columns.                                                                 |
        | dfs_month            | Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a  |
        |                      | dictionary as values. The dictionary has fields as keys (e.g. hpe, vpe) and the      |
        |                      | belonging dataframe as value with MONTHLY samples of e.g. 95th percentile and        | 
        |                      | stations as columns.                                                                 |
    """
    dsets = dset
    fields = dict()
    dfs = dict()
    dfs_day = dict()
    dfs_month = dict()  
 
    fields_def = ["east", "north", "up", "hpe", "vpe", "pos_3d", "hdop", "pdop", "vdop"]
    statistics_def = ["mean", "percentile_68", "percentile_90", "percentile_95", "rms", "std"]

    fields_cfg = config.tech[_SECTION].fields.list
    statistics_cfg = config.tech[_SECTION].statistic.list

    for field in fields_cfg:
        if field not in fields_def:
            log.fatal(f"Field {field} is not defined in 'dataset'.")
        fields.update({field: pd.DataFrame()})

    for statistic in statistics_cfg:
        if statistic not in statistics_def:
            log.fatal(f"Option '{statistic}' is not defined in 'statistic' option.")
        dfs_day.update({ statistic: copy.deepcopy(fields)})
        dfs_month.update({ statistic: copy.deepcopy(fields)})
        
    for station, dset in dsets.items():

        if dset.num_obs == 0:
            log.warn(f"Dataset '{station}' is empty.")
            continue

        # Determine topocentric coordinates (east, north, up)
        ref_pos = position.Position(
                            np.repeat(
                                  np.array([dset.meta["pos_x"], dset.meta["pos_y"], dset.meta["pos_z"]])[None, :], 
                                  dset.num_obs, 
                                  axis=0,
                            ), 
                            system="trs",
        )
        if not "enu" in dset.fields:
            dset.add_position_delta(
                name="enu",
                val=(dset.site_pos.trs - ref_pos).val,
                system="trs",
                ref_pos=ref_pos,
                write_level="operational",
            )

        # TODO: Maybe it is not necessary to introduce enu, hpe and vpe to dataset
        #      Maybe better to introduce fields in estimate stage already.
        if not "hpe" in dset.fields:
            hpe = np.sqrt(dset.enu.enu.east ** 2 + dset.enu.enu.north ** 2)
            dset.add_float("hpe", val=hpe, write_level="operational",)

        if not "vpe" in dset.fields:
            vpe = np.absolute(dset.enu.enu.up)
            dset.add_float("vpe", val=vpe, write_level="operational",)

        if not "pos_3d" in dset.fields:
            pos_3d = np.sqrt(dset.enu.enu.east ** 2 + dset.enu.enu.north ** 2 + dset.enu.enu.up ** 2)
            dset.add_float("pos_3d", val=pos_3d, write_level="operational",)


        # Determine dataframe
        df = dset.as_dataframe(fields=["enu.enu", "time.gps", "hpe", "vpe", "pos_3d", "pdop", "vdop", "hdop"])
        df = df.rename(columns={"enu_enu_0": "east", "enu_enu_1": "north", "enu_enu_2": "up"})

        if df.empty:
            continue
        else:
            # Save data in dictionaries
            dfs.update({station: df})

            for type_ in dfs_day.keys():

                df_day = _apply(df, "D", type_)
                for field in fields.keys():
                    if dfs_day[type_][field].empty:
                        dfs_day[type_][field][station] = df_day[field]
                    else:
                        dfs_day[type_][field] = pd.concat([dfs_day[type_][field], df_day[field]], axis=1)
                    dfs_day[type_][field] = dfs_day[type_][field].rename(columns={field: station})

                df_month = _apply(df, "M", type_)
                df_month.index = df_month.index.strftime("%b-%Y")
                for field in fields.keys():
                    dfs_month[type_][field][station] = df_month[field]

    return dfs, dfs_day, dfs_month


def _plot_position_error(
    dfs_day: Dict[str, pd.core.frame.DataFrame],
    dfs_month: Dict[str, pd.core.frame.DataFrame],
    figure_dir: PosixPath,
    file_vars: Dict[str, Any],
) -> None:
    """Plot horizontal and vertical position error plots 

    Args:
        dfs_day:    Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a
                    dictionary as values. The dictionary has fields as keys (e.g. hpe, vpe) and the 
                    belonging dataframe as value with DAILY samples of e.g. 95th percentile and stations as
                    columns.
        dfs_month   Dictionary with function type as keys ('mean', 'percentile_xx', 'rms', 'std') and a
                    dictionary as values. The dictionary has fields as keys (e.g. hpe, vpe) and the
                    belonging dataframe as value with MONTHLY samples of e.g. 95th percentile and stations as
                    columns.
       figure_dir:  Figure directory
    """

    ylabel_def = {
            "mean": "MEAN",
            "percentile_68": "68%",
            "percentile_90": "90%",
            "percentile_95": "95%",
            "rms": "RMS",
            "std": "STD",
    }

    options = {
        "colormap": "tab20",
        "figsize": (7, 3),
        # "grid": True,
        "marker": "o",
        "markersize": "4",
        "linestyle": "solid",
        "plot_to": "file",
        "plot_type": "plot",
        # "statistic": ["rms", "mean", "std", "min", "max", "percentile"], #TODO: Is only shown for data, which are plotted at last.
        "title": config.tech[_SECTION].title.str.upper(),
    }

    colors = (
        config.tech[_SECTION].colors.list
        if config.tech[_SECTION].colors.list
        else ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"]
    )

    # Loop over statistical solutions
    for type_ in dfs_day.keys():

        # Get used samples
        samples = dict()
        for sample in config.tech[_SECTION].samples.list:
            if "daily" == sample:
                samples["daily"] = dfs_day[type_]
            elif "monthly" == sample:
                samples["monthly"] = dfs_month[type_]
            else:
                log.fatal(f"Sample '{sample}' is not defined. Only 'daily' and/or 'monthly' can be chosen as sample.")
         
        # Loop over sampled data   
        for sample, sample_data in samples.items():

            # Loop over fields to plot
            for field in config.tech[_SECTION].fields.list:

                # Get y-range limits
                if field == "hpe":
                    ylim = config.tech[_SECTION].ylim_hpe.list
                elif field == "vpe":
                    ylim = config.tech[_SECTION].ylim_vpe.list
                elif field == "pos_3d":
                    ylim = config.tech[_SECTION].ylim_pos_3d.list
                else:
                    ylim = config.tech[_SECTION].ylim.list

                options["ylim"] = [float(ylim[0]), float(ylim[1])] if ylim else ylim

                # Generate x- and y-arrays for plotting
                x_arrays = []
                y_arrays = []
                labels = []

                for station in sample_data[field].columns:
                    #if sample == "monthly":
                    #    options.update({"xlim": "auto", "ylim": "auto"})
                    x_arrays.append(list(sample_data[field].index))
                    y_arrays.append(list(sample_data[field][station]))
                    labels.append(station.upper())

                # Generate plot
                plt = MatPlotExt()
                plt.plot(
                    x_arrays=x_arrays,
                    y_arrays=y_arrays,
                    xlabel="Time [GPS]",
                    ylabel=f"3D {ylabel_def[type_]}" if field == "pos_3d" else f"{field.upper()} {ylabel_def[type_]}",
                    y_unit="m",
                    labels=labels,
                    colors=colors,
                    figure_path=figure_dir
                    / f"plot_{type_}_{field}_{sample}_{file_vars['date']}_{file_vars['solution'].lower()}.{FIGURE_FORMAT}",
                    options=options,
                )
