"""Write report about a Galileo HAS message file analysis run

Description:
------------

asdf


"""
# Third party imports
import numpy as np
from pathlib import PosixPath
from typing import Any, Dict, List, Tuple, Union

# Midgard imports
from midgard.collections import enums
from midgard.dev import plugins
from midgard.math.unit import Unit
from midgard.plot.matplotext import MatPlotExt

# Where imports
from where import apriori
from where.lib import config
from where.lib import log
from where.lib.util import get_day_limits
from where.writers._gnss_plot import GnssPlot
from where.writers._report import Report

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])

FIGURE_FORMAT = "png"


@plugins.register
def gnss_has_report(dset: "Dataset") -> None:
    """Write report about a Galileo HAS message file analysis run

    Args:
        dset:        A dataset containing the data.
    """
    file_vars = {**dset.vars, **dset.analysis}
 
    # Generate figure directory to save figures generated for Galileo HAS message report
    figure_dir = config.files.path("output_gnss_has_report_figure", file_vars=file_vars)
    figure_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate Galileo HAS message report
    path = config.files.path("output_gnss_has_report", file_vars=file_vars)
    with config.files.open_path(path, create_dirs=True, mode="wt") as fid:
        rpt = Report(fid, rundate=dset.analysis["rundate"], path=path, description="Galileo HAS message file analysis")
        rpt.title_page()
        rpt.write_config()
        _add_to_report(dset, rpt, figure_dir)
        rpt.markdown_to_pdf()


def _add_to_report(dset: "Dataset", rpt: "Report", figure_dir: PosixPath) -> None:
    """Add figures and tables to report

    Args:
        dset:        A dataset containing the data.
        rpt:         Report object.
        figure_dir:  Figure directory.
    """
    compare_to_rinex_nav = config.tech[_SECTION].compare_to_rinex_nav.bool
    
    #
    # Galileo HAS corrections
    #
    rpt.add_text("\n# Galileo HAS corrections\n\n")
    
    for figure_path in _plot_has_correction(dset, figure_dir):
        rpt.add_figure(
                figure_path=figure_path, 
                caption=f"HAS message entry for {enums.gnss_id_to_name[figure_path.stem.split('_')[2]]}", 
                clearpage=False,
        )
        
    if "status" in dset.fields:
        for figure_path in _plot_has_clock_availability(dset, figure_dir):
            rpt.add_figure(
                    figure_path=figure_path, 
                    caption=f"HAS clock correction status for {enums.gnss_id_to_name[figure_path.stem[-1]]}", 
                    clearpage=False,
            )
          
    if "delta_radial" in dset.fields:
        for figure_path in _plot_has_orbit_availability(dset, figure_dir):
            rpt.add_figure(
                    figure_path=figure_path, 
                    caption=f"HAS orbit correction status for {enums.gnss_id_to_name[figure_path.stem[-1]]}", 
                    clearpage=True,
            )
            
    if compare_to_rinex_nav:
        
        rpt.add_text("\n# Galileo HAS corrections analysis related to broadcast navigation message\n\n")

        orbit = apriori.get(
            "orbit", 
            apriori_orbit="broadcast", 
            rundate=dset.analysis["rundate"], 
            system=tuple(dset.unique("system")),
            day_offset=0,
        )
    
        rpt.add_text("\n## GNSS signal-in-space (SIS) status\n\n")
        plt = GnssPlot(orbit.dset_edit, figure_dir, figure_format=FIGURE_FORMAT)

        # Plot GNSS SIS status (except if only Galileo system is available)
        if not (len(orbit.dset_edit.unique("system")) == 1 and orbit.dset_edit.unique("system")[0] == "E"):
            caption = "GNSS signal-in-space (SIS) status for each navigation message."
            if "E" in orbit.dset_edit.unique("system"):  # Add extra comment for Galileo
                signal, nav_type = plt.get_first_galileo_signal()
                caption += (
                    f" Galileo SIS status is given for signal '{signal.upper()}' and "
                    f"navigation message type '{nav_type}'."
                )
            rpt.add_figure(
                    figure_path=plt.plot_gnss_signal_in_space_status_overview(), 
                    caption=caption, 
                    clearpage=False,
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
                    clearpage=True,
            )
        
        rpt.add_text("\n## Galileo HAS message latency\n\n")
        
        for figure_path in _plot_latency(dset, orbit.dset_edit, figure_dir):
            if "system" in figure_path.stem:
                caption = f"Latency between HAS messages data and broadcast ephemeris for {enums.gnss_id_to_name[figure_path.stem[-1]]}"
            elif "histogram" in figure_path.stem:
                caption = f"Histogram for latency of HAS messages for {enums.gnss_id_to_name[figure_path.stem[-1]]}"
            else:
                caption = f"Latency between HAS messages data and broadcast ephemeris for satellite {figure_path.stem[-3:]}"
            rpt.add_figure(
                    figure_path=figure_path, 
                    caption=caption, 
                    clearpage=False,
            )

            
#
# PREPARE DATA
#

def _generate_latency_data(
        dset: "Dataset", 
        orbit: "Dataset",
) -> Tuple[Dict[str, List[np.ndarray]], Dict[str, List[np.ndarray]], List[str], List[Tuple[Any]]]:
    """Generate latency data based on time difference between HAS and broadcast navigation messages

    Args:
        dset:     A dataset containing the HAS message data.
        orbit:    A dataset continaing the broadcast navigation message data
        
    Returns:
        Tuple with x-arrays and y-arrays dictionaries with satellite keys and HAS reception time and latency
        as values
        
       | Element             | Type                         | Description                                                     |
       |---------------------|------------------------------|-----------------------------------------------------------------|
       | x_arrays            | Dict[str, List[np.ndarray]]  | Time entries for each satellite used as x-array                 |
       | y_arrays            | Dict[str, List[np.ndarray]]  | Latency entries for each satellites used as y-array             |
       | labels              | Dict[str, List[np.ndarray]]  | Satellite labels                                                |
       | missing_satellites  | List[str]                    | Satellites with no HAS messages in relation to broadcast        |
       |                     |                              | navigation message                                              |
       | missing_messages    | List[Tuple[Any]]             | Missing HAS messages in relation to broadcast navigation        |
       |                     |                              | messages given as tuple (SAT, IODE) or (SAT, IODE, TIME)        |
    """
       
    # Generate empty data dictionaries
    missing_satellites = list()
    missing_messages = list()
    x_arrays = dict()
    y_arrays = dict()
    labels = dict()
    
    for sat in orbit.unique("satellite"):
        
        if not sat in dset.unique("satellite"):
            if sat[0] in dset.unique("system"):
                log.warn(f"Satellite {sat} in broadcast navigation message is not given in HAS messages.")
                missing_satellites.append(sat)
            continue
        
        x_arrays.update({sat: list()})
        y_arrays.update({sat: list()})
        labels.update({sat: list()})  
        
    # Fill data dictionaries for each satellite
    for sat in x_arrays.keys():
              
        idx_sat = orbit.filter(satellite=sat)
        
        # Generate data to plot
        for iode in set(orbit.iode[idx_sat]):
        
            idx_nav = orbit.filter(satellite=sat, iode=iode)
            idx_has = dset.filter(satellite=sat, gnssiod=iode)
            labels[sat].append(iode)
            
            if not dset.time[idx_has].size == 0:

                x_tmp = np.array([])                    
                y_tmp = np.array([])
                
                # Loop over time is necessary due to that an IOD can be used several times on a day
                for time in orbit.transmission_time[idx_nav]:
                                        
                    diff = (dset.time.gps.mjd[idx_has] - time.gps.mjd) * Unit.day2minute
                    idx = abs(diff) < 1166.667 # (1166.667 min ~= 70000 s) Skip times from 2nd IOD -> TODO: Better solution to handle that?
                
                    if dset.time.gps.datetime[idx_has][idx].size == 0:
                        log.debug(f"For satellite {sat}, IOD {iode} and broadcast transmission time {time.gps.isot} does not exist HAS messages.")
                        missing_messages.append({sat, iode, time})
                        continue
                    
                    x_tmp = np.hstack([x_tmp, dset.time.gps.datetime[idx_has][idx]])
                    y_tmp = np.hstack([y_tmp, diff[idx]])
                                                   
                    log.debug(f"Time difference between HAS and broadcast navigation message wit IOD {iode:3.0f} is "
                             f"between {min(diff):7.0f} min and {max(diff):7.0f} min for satellite {sat}")
                
                if x_tmp.size == 0: #TODO: check in which cases this happens?
                    continue
                
                x_arrays[sat].append(x_tmp)   
                y_arrays[sat].append(y_tmp)
                   
            else:
                log.debug(f"For satellite {sat} and IOD {iode} does not exist HAS messages.")
                missing_messages.append({sat, iode})
                                
    return x_arrays, y_arrays, labels, missing_satellites, missing_messages


#
# PLOT FUNCTIONS
#

def _get_galileo_has_clock_status_data(
            dset: "Dataset",
            status: int, 
            system: Union[str, None]=None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Get time and satellite data based on a given Galileo HAS clock correction status 

    The Galileo HAS clock correction status can be:

     | CODE | CLOCK STATUS                 | 
     |------|------------------------------|
     |   0  | data ok                      |
     |   1  | data not available           | 
     |   2  | satellite shall not be used  |

    Args:
       dset:          A dataset containing the data.
       status:        Galileo HAS clock correction status 
       system:        Get only Galileo HAS clock correction status for selected system identifier.

    Returns:
        Tuple with time and satellite data for a given Galileo HAS clock correction status and ordered by satellite
    """
    # TODO: order of satellites is not correct. Maybe to save the data in a dataframe could help?

    # Generate x- and y-axis data
    time = []
    satellite = []

    for sat in sorted(dset.unique("satellite"), reverse=True):
        
        if system: # Skip GNSS, which are not selected
            if not sat.startswith(system):
                continue
            
        idx = dset.filter(satellite=sat)
        idx_status = dset.status[idx] == status

        time.extend(dset.time.gps.datetime[idx][idx_status])
        satellite.extend(dset.satellite[idx][idx_status])

    return time, satellite


def _get_galileo_has_orbit_status_data(
                dset: "Dataset",
                data_available: bool, 
                system: str,
    ) -> Tuple[np.ndarray, np.ndarray]:
    """Get time and satellite data based on a given Galileo HAS orbit correction status 

    The Galileo HAS orbit correction status can be accessed via Not a Number (NaN). If data are NaN, then orbit
    correction are not available otherwise they are available.

     | Type  | ORBIT STATUS                 | 
     |-------|------------------------------|
     | float | data available               |
     | NaN   | data not available           | 
     
    TODO: Has to be tested.

    Args:
       dset:            A dataset containing the data.
       data_available:   Depending on this flag it is checked, if either orbit correction are available or not
       system:           Get only Galileo HAS clock orbit status for selected system identifier.

    Returns:
        Tuple with time and satellite data for a given Galileo HAS orbit correction status and ordered by satellite
    """
    
    time = []
    satellite = []

    for sat in sorted(dset.unique("satellite"), reverse=True):
        
        if system: # Skip GNSS, which are not selected
            if not sat.startswith(system):
                continue
            
        idx = dset.filter(satellite=sat)
        idx_status = np.isnan(dset.delta_radial[idx])
        idx_status = ~idx_status if data_available else idx_status

        time.extend(dset.time.gps.datetime[idx][idx_status])
        satellite.extend(dset.satellite[idx][idx_status])
  
    return time, satellite


def _plot_has_clock_availability(
                        dset: "Dataset",
                        figure_dir:  PosixPath,
                        figure_name: str="plot_galileo_has_clock_availability_{system}.{FIGURE_FORMAT}",

) -> List[PosixPath]:
    """Plot Galileo HAS clock correction availability
    
    Args:
        dset:        A dataset containing the data.
        figure_dir:  Figure directory.
        figure_name: File name of figure.
    
    Returns:
        List with figure path for Galileo HAS clock correction availability plots. File name ends with GNSS 
        identifier (e.g. 'E', 'G') and observation type, for example 'plot_galileo_has_clock_availability_G.png'.
            
    """
    figure_paths = list()
    colors = ["green", "yellow", "red"]
    labels = ["data ok", "data not available", "data shall not be used"]
    status_def = [0, 1, 2]
    
    for system in dset.unique("system"):
        
        figure_path = figure_dir / figure_name.replace("{FIGURE_FORMAT}", FIGURE_FORMAT).replace("{system}", system)
        
        # Generate time and satellite data for given Galileo HAS clock correction status
        x_arrays = []
        y_arrays = []
        for status in status_def:   
            time, satellite = _get_galileo_has_clock_status_data(dset, status, system)
            x_arrays.append(time)
            y_arrays.append(satellite)
    
        # Limit x-axis range to rundate
        day_start, day_end = get_day_limits(dset)
    
        # Generate plot
        plt = MatPlotExt()
        plt.plot(
            x_arrays=x_arrays,
            y_arrays=y_arrays,
            xlabel="Time [GPS]",
            ylabel="Satellite",
            y_unit="",
            labels=labels,
            colors=colors,
            figure_path=figure_path,
            options={
                "figsize": (7, 8),
                "marker": "s",
                "marksersize": 10,
                "legend_ncol": 4,
                "legend_location": "bottom",
                "plot_to": "file",
                "plot_type": "scatter",
                "tick_labelsize": ("y", 8),  # specify labelsize for y-axis
                "title": "HAS clock correction availability",
                "xlim": [day_start, day_end],
            },
        )
        # TODO: Legend has to be improved. Old configuration:
        # figsize = (7, 10)
        # loc="lower center",
        # bbox_to_anchor=(0.5, -0.01),
        # frameon=True,
        
        # Add figure path
        figure_paths.append(figure_path)
  
    return figure_paths


def _plot_has_orbit_availability(
                        dset: "Dataset", 
                        figure_dir: PosixPath,
                        figure_name: str="plot_galileo_has_orbit_availability_{system}.{FIGURE_FORMAT}",

) -> List[PosixPath]:
    """Plot Galileo HAS orbit correction availability
    
    Args:
        dset:        A dataset containing the data.
        figure_dir:  Figure directory.
        figure_name: File name of figure.
    
    Returns:
        List with figure path for Galileo HAS orbit correction availability plots. File name ends with GNSS 
        identifier (e.g. 'E', 'G') and observation type, for example 'plot_galileo_has_orbit_availability_G.png'.
            
    """
    figure_paths = list()
    colors = ["green", "red"]
    labels = ["data ok", "data not available"]
    
    for system in dset.unique("system"):
                  
        figure_path = figure_dir / figure_name.replace("{FIGURE_FORMAT}", FIGURE_FORMAT).replace("{system}", system)
        
        # Generate time and satellite data for given Galileo HAS orbit correction status
        x_arrays = []
        y_arrays = []
        
        for data_available in [True, False]:   
            time, satellite = _get_galileo_has_orbit_status_data(dset, data_available, system)
            x_arrays.append(time)
            y_arrays.append(satellite)
    
        # Limit x-axis range to rundate
        day_start, day_end = get_day_limits(dset)
    
        # Generate plot
        plt = MatPlotExt()
        plt.plot(
            x_arrays=x_arrays,
            y_arrays=y_arrays,
            xlabel="Time [GPS]",
            ylabel="Satellite",
            y_unit="",
            labels=labels,
            colors=colors,
            figure_path=figure_path,
            options={
                "figsize": (7, 8),
                "marker": "s",
                "marksersize": 10,
                "legend_ncol": 4,
                "legend_location": "bottom",
                "plot_to": "file",
                "plot_type": "scatter",
                "tick_labelsize": ("y", 8),  # specify labelsize for y-axis
                "title": "HAS orbit availability",
                "xlim": [day_start, day_end],
            },
        )
        # TODO: Legend has to be improved. Old configuration:
        # figsize = (7, 10)
        # loc="lower center",
        # bbox_to_anchor=(0.5, -0.01),
        # frameon=True,
        
        # Add figure path
        figure_paths.append(figure_path)
  
    return figure_paths
  

def _plot_has_correction(
                        dset: "Dataset",
                        figure_dir: PosixPath,
                        figure_name: str="plot_{solution}.{FIGURE_FORMAT}",

) -> List[PosixPath]:
    """Plot Galileo HAS correction results
    
    Args:
        dset:        A dataset containing the data.
        figure_dir:  Figure directory.
        figure_name: File name of figure.
    
    Returns:
        List with figure path for Galileo HAS correction plots. File name includes GNSS identifier
        (e.g. 'E', 'G') and field name, for example 'plot_field_G_age_of_ephemeris.png'.
    """          
    figure_paths = list()
    
    gnss_plt = GnssPlot(dset, figure_dir, figure_format=FIGURE_FORMAT)
           
    # Define fields to plot
    fields_without_labels={"age_of_data", "iod"}
    remove_fields = {"satellite", "signal", "status", "system", "time", "tom"}
    fields = set(dset.fields) - remove_fields
    
    # Get file name suffix
    if "delta_radial" in dset.fields:
        suffix = "_orb"
    elif "multiplier" in dset.fields:
        suffix = "_clk"
    elif "code_bias" in dset.fields:
        suffix = "_cb"
    elif "phase_bias" in dset.fields:
        suffix = "_pb"
    
    #
    # Plot HAS message fields 
    #         
    for field in fields:        
        if field in fields_without_labels:
            figure_paths = figure_paths + gnss_plt.plot_field(field, use_labels=False, suffix=suffix) # without satellite labels
        else:
            figure_paths = figure_paths + gnss_plt.plot_field(field, suffix=suffix) # with satellite labels
            
    return figure_paths



def _plot_latency(dset: "Dataset", orbit: "Dataset", figure_dir: PosixPath)  -> List[PosixPath]:
    """Generate latency data based on time difference between HAS and broadcast navigation messages

    Args:
        dset:        A dataset containing the data.
        figure_dir:  Figure directory.

    Returns:
        List with figure path. File name ends with GNSS satellite identifier (e.g. 'E01', 'G02'), for example 
        'plot_latency_E01.png' or GNSS system identifier (e.g. E, G), for example 'plot_latency_system_E.png'.
    """
    x_arrays, y_arrays, labels, missing_satellites, missing_messages = _generate_latency_data(dset, orbit)
    figure_paths = list()
    file_suffix = dset.vars["file_key"].split("_")[-1]
    
    # Generate dataset with complete satellite data
    for sys in dset.unique("system"):
        x_arrays_sys = list()
        y_arrays_sys = list()
        labels_sys = list()
        
        for sat in sorted(x_arrays.keys()):
            if not sat[0] == sys: # Skip satellite not belonging to chosen GNSS
                continue
            labels_sys.append(sat)
            x_tmp = list()
            y_tmp = list()
            
            for x_array, y_array in zip(x_arrays[sat], y_arrays[sat]):
                x_tmp.extend(list(x_array))
                y_tmp.extend(list(y_array))
                
            x_arrays_sys.append(x_tmp)
            y_arrays_sys.append(y_tmp)
        
        # Plot latency for all satellites
        figure_path = figure_dir / f"plot_latency_system_{file_suffix}_{sys}.{FIGURE_FORMAT}"
        figure_paths.append(figure_path)
             
        plt = MatPlotExt()
        plt.plot(
            x_arrays=x_arrays_sys,
            y_arrays=y_arrays_sys,
            xlabel="Reception time of HAS messages",
            ylabel="Latency to navigation message",
            y_unit="minute",
            labels=labels_sys,
            figure_path=figure_path,
            options={
                "colormap": "jet",
                "figsize": (7, 7),
                "legend": True,
                "legend_ncol": 6,
                "legend_location": "bottom",
                "plot_to": "file",
                "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
            },
        )
        
        # Plot histogram
        data = list()
        for y_array in y_arrays_sys:
            data.extend(y_array)
            
        figure_paths.append(_plot_histogram_latency(data, figure_dir, sys, file_suffix))
  
    # Plot latency for each satellite
    for sat in sorted(x_arrays.keys()):
        
        figure_path = figure_dir / f"plot_latency_{file_suffix}_{sat}.{FIGURE_FORMAT}"
        figure_paths.append(figure_path)
        
        plt = MatPlotExt()
        plt.plot(
            x_arrays=x_arrays[sat],
            y_arrays=y_arrays[sat],
            xlabel="Reception time of HAS messages",
            ylabel="Latency to navigation message",
            y_unit="minute",
            labels=labels[sat],
            figure_path=figure_path,
            options={
                "colorbar": True,
                "colorbar_label": "IOD",
                "colormap": "jet",
                "figsize": (7, 4),
                "plot_to": "file",
                "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
                "title": f"Satellite {sat}",
            },
        )
                         
    return figure_paths


#
# HISTOGRAM PLOTS
#
def _plot_histogram_subplot(data: List[float], axis: "AxesSubplot", system: str):
    """Plot histogram subplots

    Args:
       data:    Data to plot.
       axis:    Subplot axes.
       system:  GNSS system identifier (e.g. E, G, ...)
    """
    axis.hist(data, density=True, bins=50)
    axis.set(xlabel="Latency [min]", ylabel="Frequency")
    axis.set_title(f"{enums.gnss_id_to_name[system]}")
    mean = np.mean(data)
    std = np.std(data)
    axis.text(
        0.98,
        0.98,
        f"$mean={mean:5.3f} \\pm {std:5.3f}$\n#data={len(data)}",
        horizontalalignment="right",
        verticalalignment="top",
        transform=axis.transAxes,
    )


def _plot_histogram_latency(data: List[float], figure_dir: PosixPath, sys: str, file_suffix: str) -> PosixPath:
    """Generate latency histogram plot

    Args:
        data:        Data for histogram plot
        figure_dir:  Figure directory.
        sys:         GNSS system
        file_suffix: File suffix, which defines Galileo HAS correction type (e.g. 'clk', 'orb')

    Returns:
        Figure path.
    """
    import matplotlib.pyplot as plt
    figure_path = figure_dir / f"plot_histogram_latency_{file_suffix}_{sys}.{FIGURE_FORMAT}"

    fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, squeeze=False)
    _plot_histogram_subplot(data, axes.flatten()[0], sys)
    plt.savefig(figure_path, dpi=200)
    plt.clf()  # clear the current figure

    return figure_path          
            
        

    

