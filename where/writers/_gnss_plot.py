"""Class for GNSS plots

Description:
------------
TODO

"""
# Standard liberay imports
from dataclasses import dataclass
from collections import namedtuple
import pathlib
from typing import List, Tuple, Union

# External liberary imports
import numpy as np

# Midgard imports
from midgard.collections import enums
from midgard.dev import exceptions
from midgard.gnss import gnss
from midgard.math.unit import Unit
from midgard.plot.matplotext import MatPlotExt

# Where imports
from where.data import dataset3 as dataset
from where.lib import config
from where.lib import log
from where.lib.util import get_day_limits
    


@dataclass
class PlotField:
    """A convenience class for defining a necessary plotting parameters

    Args:
        ylabel:   y-axis label description
        unit:     Unit used for plotting of y-axis
        figsize:  Figure size
    """
    ylabel: str
    unit: str
    figsize: Tuple[float, float]
    

class GnssPlot:
    """Class for GNSS plots
    """

    def __init__(
        self, 
        dset: "Dataset",
        figure_dir: pathlib.PosixPath,
        figure_format: str="png",
    ) -> None:
        """Set up a new GNSS plot object

        Args:
            dset:          A dataset containing the data.
            figure_dir:    Figure directory.
            figure_format: Figure format.
        """
        self.dset = dset
        self.figure_dir = figure_dir
        self.figure_format = figure_format
        
        
    def plot_dop(self,
                  figure_name: str="plot_dop.{FIGURE_FORMAT}",
    ) -> pathlib.PosixPath:
        """Plot DOP
    
        Args:
            figure_name: File name of figure.
        """
        figure_path = self.figure_dir / figure_name.replace("{FIGURE_FORMAT}", self.figure_format)
        log.debug(f"Plot {figure_path}.")
    
        plt = MatPlotExt()
        plt.plot(
            x_arrays=[
                self.dset.time.gps.datetime,
                self.dset.time.gps.datetime,
                self.dset.time.gps.datetime,
                self.dset.time.gps.datetime,
                self.dset.time.gps.datetime,
            ],
            y_arrays=[self.dset.gdop, self.dset.pdop, self.dset.vdop, self.dset.hdop, self.dset.tdop],
            xlabel="Time [GPS]",
            ylabel="Dilution of precision",
            y_unit="",
            labels=["GDOP", "PDOP", "VDOP", "HDOP", "TDOP"],
            figure_path=figure_path,
            options={
                "figsize": (7, 4), 
                "legend": True, 
                "plot_to": "file",
            },
        )
        
        return figure_path
        

    def plot_epoch_by_epoch_difference(
                            self, 
                            figure_name: str="plot_epoch_by_epoch_difference_{solution}.{FIGURE_FORMAT}",

    ) -> List[pathlib.PosixPath]:
        """Plot epoch by epoch difference of observations
        
        Args:
            figure_name: File name of figure.
        
        Returns:
            List with figure path for linear combinations depending on GNSS. File name ends with GNSS identifier
            (e.g. 'E', 'G') and observation type, for example 'plot_geometry_free_code_G_C1C_C2X.png'.
        """
        figure_paths = list()
        
        for field in self.dset.diff_epo.fields:

            for sys in sorted(self.dset.meta["obstypes"].keys()):

                x_arrays = []
                y_arrays = []
                labels = []  
                
                # Skip if nothing to plot
                idx_sys = self.dset.filter(system=sys)
                if np.all(np.isnan(self.dset.diff_epo[field][idx_sys])):
                    continue
                                  
                figure_path = self.figure_dir / figure_name.replace("{solution}", f"{sys}_{field}").replace("{FIGURE_FORMAT}", self.figure_format)
                figure_paths.append(figure_path)
                log.debug(f"Plot {figure_path}.")
    
                for sat in sorted(self.dset.unique("satellite")):
                    if not sat.startswith(sys):
                        continue
                    idx = self.dset.filter(satellite= sat)
                    x_arrays.append(self.dset.time.gps.datetime[idx])
                    y_arrays.append(self.dset.diff_epo[field][idx])
                    labels.append(sat)  
                
                # Plot scatter plot
                plt = MatPlotExt()
                plt.plot(
                    x_arrays=x_arrays,
                    y_arrays=y_arrays,
                    xlabel="Time [GPS]",
                    ylabel=f"Epoch by epoch difference ({field})",
                    y_unit="m",
                    labels=labels,
                    figure_path=figure_path,
                    options={
                        "figsize": (7, 6),
                        "legend": True,
                        "legend_ncol": 6,
                        "legend_location": "bottom",
                        "plot_to": "file", 
                        "plot_type": "scatter",
                        "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
                    },
                )
                
        return figure_paths
    
    
    def plot_field(
                    self, 
                    field: str,
                    collection: Union[str, None] = None, 
                    figure_name: str="plot_field_{solution}.{FIGURE_FORMAT}",
                    use_labels: bool=True,
                    suffix: str="",

    ) -> List[pathlib.PosixPath]:
        """Plot field
        
        Args:
            collection:  Collection name.
            field:       Field name.
            figure_name: File name of figure.
            use_labels:  Use labels together with legend.
            suffix:      File name suffix
        
        Returns:
            List with figure path for depending on GNSS. File name ends with GNSS identifier (e.g. 'E', 'G') and field 
            name, for example 'plot_field_G_C1C.png'.
        """
        figure_paths = list()
        fieldname = f"{collection}.{field}" if collection else field
        
        plot_def = {
                "age_of_ephemeris": PlotField("Age of ephemeris", "s", (7, 6)),
                "bgd_e1_e5a": PlotField("BGD(E1,E5a)", "ns", (7, 6)),
                "bgd_e1_e5a_diff": PlotField("BGD(E1,E5a) - DCB(C1C,C5Q)", "ns", (7, 6)),
                "bgd_e1_e5a_diff_mean": PlotField("BGD(E1,E5a) - DCB(C1C,C5Q)", "ns", (7, 6)),
                "bgd_e1_e5b": PlotField("BGD(E1,E5b)", "ns", (7, 6)),
                "bgd_e1_e5b_diff": PlotField("BGD(E1,E5b) - DCB(C1C,C7Q)", "ns", (7, 6)),
                "bgd_e1_e5b_diff_mean": PlotField("BGD(E1,E5b) - DCB(C1C,C7Q)", "ns", (7, 6)),
                "clk_diff_with_dt_mean": PlotField("Clock correction difference $\Delta t$ (mean)", "m", (7, 6)),
                "delay.gnss_earth_rotation_drift": PlotField("Earth rotation drift", "m/s", (7, 6)),
                "delay.gnss_ionosphere": PlotField("Ionospheric delay", "m", (7, 6)),
                "delay.gnss_range": PlotField("Range", "m", (7, 6)),
                "delay.gnss_range_rate": PlotField("Range rate", "m/s", (7, 6)),
                "delay.gnss_relativistic_clock": PlotField("Relativistic clock", "m", (7, 6)),
                "delay.gnss_relativistic_clock_rate": PlotField("Relativistic clock rate", "m/s", (7, 6)),
                "delay.gnss_satellite_clock": PlotField("Satellite clock", "m", (7, 6)),
                "delay.gnss_satellite_clock_rate": PlotField("Satellite clock rate", "m/s", (7, 6)),
                "delay.gnss_total_group_delay": PlotField("Total group delay", "m", (7, 6)),
                "delay.troposphere_radio": PlotField("Troposphere delay", "m", (7, 6)),
                "orb_diff_3d": PlotField("3D orbit error", "m", (7, 6)),
                "sisre": PlotField("SISE", "m", (7, 6)),
                "sisre_orb": PlotField("orbit-only SISE", "m", (7, 6)),
                "tgd": PlotField("TGD(L1,L2)", "ns", (7, 6)),
                "tgd_diff": PlotField("TGD(L1,L2) - DCB(C1W,C2W)", "ns", (7, 6)),
                "tgd_diff_mean": PlotField("TGD(L1,L2) - DCB(C1W,C2W)", "ns", (7, 6)),
                "tgd_b1_b3": PlotField("TGD(B1,B3)", "ns", (7, 6)),
                "tgd_b1_b3_diff": PlotField("TGD(B1,B3) - DCB(C2I,C6I)", "ns", (7, 6)),
                "tgd_b1_b3_diff_mean": PlotField("TGD(B1,B3) - DCB(C2I,C6I)", "ns", (7, 6)),
                "tgd_b2_b3": PlotField("TGD(B2,B3)", "ns", (7, 6)),
                "tgd_b2_b3_diff": PlotField("TGD(B2,B3) - DCB(C7I,C6I)", "ns", (7, 6)),
                "tgd_b2_b3_diff_mean": PlotField("TGD(B2,B3) - DCB(C7I,C6I)", "ns", (7, 6)),

                # Galileo HAS plots
                "age_of_data": PlotField("Age of data", "s", (7, 5)),
                "code_bias": PlotField("HAS code bias", "m", (7, 7)),
                "delta_clock_c0": PlotField("Δclock-c0", "m", (7, 7)),
                "delta_cross_track": PlotField("Δcross-track", "m", (7, 7)),
                "delta_in_track": PlotField("Δin-track", "m", (7, 7)),
                "delta_radial": PlotField("Δradial", "m", (7, 7)),
                "gnssid": PlotField("HAS GNSS ID", "", (7, 7)),
                "gnssiod": PlotField("HAS GNSS IOD", "", (7, 7)),
                "iod": PlotField("HAS IOD", "", (7, 5)),
                "multiplier": PlotField("Multiplier for Δclock-c0", "", (7, 7)),
                "phase_bias": PlotField("HAS phase bias", "cycle", (7, 7)),
                "validity": PlotField("HAS validity interval", "s", (7, 7)),

        }

        plot_def_qzss = {
                "tgd_diff": PlotField("TGD(L1,L2) - DCB(C1X,C2X)", "ns", (7, 6)),
                "tgd_diff_mean": PlotField("TGD(L1,L2) - DCB(C1X,C2X)", "ns", (7, 6)),
        }


        #+Change ylim
        ylim = []
        #if field == "clk_diff_with_dt_mean":
        #    ylim = [-0.8, 1.3]
        #elif field == "sisre":
        #    ylim = [0.0, 1.3]
        #elif field == "sisre_orb":
        #    ylim = [0.0, 0.5]
        #-Change ylim

        for sys in sorted(self.dset.unique("system")):

            x_arrays = []
            y_arrays = []
            labels = []  
            
            # Skip if nothing to plot
            idx_sys = self.dset.filter(system=sys)
            if np.all(np.isnan(self.dset[fieldname][idx_sys])):
                continue
                              
            figure_path = self.figure_dir / figure_name.replace("{solution}", f"{sys}_{field}{suffix}").replace("{FIGURE_FORMAT}", self.figure_format)
            figure_paths.append(figure_path)
            log.debug(f"Plot {figure_path}.")

            for sat in sorted(self.dset.unique("satellite")):
                if not sat.startswith(sys):
                    continue
                idx = self.dset.filter(satellite= sat)
                labels.append(sat)  
                x_arrays.append(self.dset.time.gps.datetime[idx])
                
                # Convert y_array to defined unit
                if fieldname in plot_def.keys() and plot_def[fieldname].unit:         
                    y_array = self._convert_to_unit(self.dset[fieldname][idx], self.dset.unit(fieldname)[0], plot_def[fieldname].unit) 
                else:
                    y_array = self.dset[fieldname][idx]
                
                y_array = y_array
                y_arrays.append(y_array) 
 
            if fieldname in plot_def.keys():
                if sys == "J" and fieldname in ["tgd_diff", "tgd_diff_mean"]:
                    ylabel = plot_def_qzss[fieldname].ylabel
                else:
                    ylabel = plot_def[fieldname].ylabel
            else:
                ylabel = f"Field ({field})"

            # Plot scatter plot
            plt = MatPlotExt()
            plt.plot(
                x_arrays=x_arrays,
                y_arrays=y_arrays,
                xlabel="Time [GPS]",
                ylabel=ylabel,
                y_unit=plot_def[fieldname].unit if fieldname in plot_def.keys() else "",
                labels=labels if use_labels else None,
                figure_path=figure_path,
                options={
                    "figsize": plot_def[fieldname].figsize if fieldname in plot_def.keys() else (7,6),
                    "legend": True,
                    "legend_ncol": 6,
                    "legend_location": "bottom",
                    "plot_to": "file", 
                    "plot_type": "scatter",
                    #"statistic": ["rms", "mean", "std", "min", "max", "percentile"],  #TODO: Only for last satellite solution the statistic is plotted.
                    "xlim": [min(self.dset.time.gps.datetime), max(self.dset.time.gps.datetime)],
                    "ylim": ylim,
                },
            )
                
        return figure_paths
    
    
    def plot_gnss_signal_in_space_status(
                        self,
                        figure_name: str="plot_signal_in_space_status_{solution}.{FIGURE_FORMAT}",
    ) -> List[pathlib.PosixPath]:
        """Generate GNSS Signal-in-Space (SIS) status plot for each GNSS based on GNSS signal health status (SHS) 
        and for Galileo in addition based on SIS Accuracy (SISA) and data validity status (DVS) given in RINEX 
        navigation file.
    
        The SIS status can be:
    
         | CODE | SIS STATUS      | PLOTTED COLOR | DESCRIPTION                     |
         |------|-----------------|---------------|---------------------------------|
         |   0  | healthy         |         green | SIS status used by all GNSS     |
         |   1  | marginal (SISA) |        yellow | SIS status only used by Galileo |
         |   2  | marignal (DVS)  |        orange | SIS status only used by Galileo |
         |   3  | unhealthy       |           red | SIS status used by all GNSS     |
    
        Args:
            figure_name: File name of figure.
        """
        figure_paths = list()
        colors = {
                "C": ["green", "red"],
                "E": ["green", "yellow", "orange", "red"],
                "G": ["green", "red"],
                "I": ["green", "red"],
                "J": ["green", "red"],
        }
        labels = {
                "C": ["healthy", "unhealthy"],
                "E": ["healthy", "marginal (sisa)", "marginal (dvs)", "unhealthy"],
                "G": ["healthy", "unhealthy"],
                "I": ["healthy", "unhealthy"],
                "J": ["healthy", "unhealthy"],
        }
        status_def = {
                "C": [0, 3],
                "E": [0, 1, 2, 3],
                "G": [0, 3],
                "I": [0, 3],
                "J": [0, 3],
        }
        
        # Loop over GNSSs
        for sys in self.dset.unique("system"):

            signals = self._select_gnss_signal(sys)
            
            # Generate plot for each given navigation type and in case of Galileo in addition for each signal 
            # (e.g. E1, E5a, E5b)
            for signal, nav_type in sorted(signals.items()):
                x_arrays = []
                y_arrays = []
                
                solution = f"{enums.gnss_id_to_name[sys].value.lower()}_{signal}" if signal else f"{enums.gnss_id_to_name[sys].value.lower()}"
                figure_path=self.figure_dir / figure_name.replace("{solution}", solution).replace("{FIGURE_FORMAT}", self.figure_format)
                figure_paths.append(figure_path)
                
                for status in status_def[sys]:
                    if sys == "E":
                        time, satellite = self._get_gnss_signal_in_space_status_data(status, signal, system=sys)
                    else:
                        time, satellite = self._get_gnss_signal_in_space_status_data(status, system=sys)
        
                    x_arrays.append(time)
                    y_arrays.append(satellite)
        
                # Limit x-axis range to rundate
                day_start, day_end = get_day_limits(self.dset)
                
                # Title
                if sys == "E":
                    title = f"{enums.gnss_id_to_name[sys].value.upper()} signal-in-space status for signal {signal.upper()} ({nav_type})"
                else:    
                    title = f"{enums.gnss_id_to_name[sys].value.upper()} signal-in-space status ({nav_type})"
                    
                # Generate plot
                plt = MatPlotExt()
                plt.plot(
                    x_arrays=x_arrays,
                    y_arrays=y_arrays,
                    xlabel="Time [GPS]",
                    ylabel="Satellite",
                    y_unit="",
                    labels=labels[sys],
                    colors=colors[sys],
                    figure_path=figure_path,
                    options={
                        "figsize": (7, 8),
                        "marker": "s",
                        "marksersize": 10,
                        "legend_ncol": 4,
                        "legend_location": "bottom",
                        "plot_to": "file",
                        "plot_type": "scatter",
                        "title": title,
                        "xlim": [day_start, day_end],
                    },
                )
            
        return figure_paths
    
    
    def plot_gnss_signal_in_space_status_overview(
                    self,
                    figure_name: str="plot_gnss_signal_in_space_status.{FIGURE_FORMAT}",
    ) -> pathlib.PosixPath:
        """Generate GNSS Signal-in-Space (SIS) status overview plot based on SIS status given in RINEX navigation file
    
        The SIS status can be:
    
         | CODE | SIS STATUS      | PLOTTED COLOR | DESCRIPTION                     |
         |------|-----------------|---------------|---------------------------------|
         |   0  | healthy         |         green | SIS status used by all GNSS     |
         |   1  | marginal (SISA) |        yellow | SIS status only used by Galileo |
         |   2  | marignal (DVS)  |        orange | SIS status only used by Galileo |
         |   3  | unhealthy       |           red | SIS status used by all GNSS     |
    
        Args:
            figure_name: File name of figure.
        """
        figure_path = self.figure_dir / figure_name.replace("{FIGURE_FORMAT}", self.figure_format)
        colors = ["green", "yellow", "orange", "red"]
        labels = ["healthy", "marginal (sisa)", "marginal (dvs)", "unhealthy"]
        status_def = [0, 1, 2, 3]
        signal = None
    
        # Select only one Galileo signal
        # Note: Navigation message for Galileo can include I/NAV and F/NAV messages for different signals (E1, E5a, E5b).
        #       For plotting we choose only one of them.
        if "E" in self.dset.unique("system"):
            signal, _ = self.get_first_galileo_signal()
    
        # Generate time and satellite data for given SIS status
        x_arrays = []
        y_arrays = []
        for status in status_def:
            time, satellite = self._get_gnss_signal_in_space_status_data(status, signal)
            x_arrays.append(time)
            y_arrays.append(satellite)
    
        # Limit x-axis range to rundate
        day_start, day_end = get_day_limits(self.dset)
    
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
                "figsize": (7, 12),
                "marker": "s",
                "marksersize": 10,
                "legend_ncol": 4,
                "legend_location": "bottom",
                "plot_to": "file",
                "plot_type": "scatter",
                "tick_labelsize": ("y", 7),  # specify labelsize 7 for y-axis
                "title": "GNSS signal-in-space status",
                "xlim": [day_start, day_end],
            },
        )
        # TODO: Legend has to be improved. Old configuration:
        # figsize = (7, 10)
        # loc="lower center",
        # bbox_to_anchor=(0.5, -0.01),
        # frameon=True,
        
        return figure_path
    
             
    def plot_linear_combinations(
                            self, 
                            figure_name: str="plot_{solution}.{FIGURE_FORMAT}",

    ) -> List[pathlib.PosixPath]:
        """Plot linear combinations of observations
        
        Args:
            figure_name: File name of figure.
        
        Returns:
            List with figure path for linear combinations depending on GNSS. File name ends with GNSS identifier
            (e.g. 'E', 'G') and observation type, for example 'plot_geometry_free_code_G_C1C_C2X.png'.
        """
        FigureInfo = namedtuple("FigureInfo", ["file_path", "name", "system", "obstype"])
        
        figure_info = list()
    
        for field in self.dset.lin.fields:

            # Melbourne-Wübbena linear combination
            if field == "melbourne_wuebbena":
                 
                for sys in sorted(self.dset.meta["obstypes"].keys()):
                    x_arrays = []
                    y_arrays = []
                    labels = []                    
                    solution = f"{field}_{sys}_{'_'.join(self.dset.meta['linear_combination'][field][sys])}"
                    figure_path = self.figure_dir / figure_name.replace("{solution}", solution).replace("{FIGURE_FORMAT}", self.figure_format)
                    log.debug(f"Plot {figure_path}.")
                    figure_info.append(FigureInfo(
                                        figure_path, 
                                        field, 
                                        sys, 
                                        self.dset.meta['linear_combination'][field][sys]),
                    )
        
                    for sat in sorted(self.dset.unique("satellite")):
                        if not sat.startswith(sys):
                            continue
                        idx = self.dset.filter(satellite= sat)
                        x_arrays.append(self.dset.time.gps.datetime[idx])
                        y_arrays.append(self.dset.lin[field][idx])
                        labels.append(sat)  
                    
                    # Plot scatter plot
                    plt = MatPlotExt()
                    plt.plot(
                        x_arrays=x_arrays,
                        y_arrays=y_arrays,
                        xlabel="Time [GPS]",
                        ylabel="Melbourne-Wübbena",
                        y_unit="m",
                        labels=labels,
                        figure_path=figure_path,
                        options={
                            "figsize": (7, 6),
                            "legend": True,
                            "legend_ncol": 6,
                            "legend_location": "bottom",
                            "plot_to": "file", 
                            "plot_type": "scatter",
                            "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
                            "xlim": "fit_to_data",
                        },
                    )
                
            # Code-multipath linear combination
            elif field == "code_multipath_f1" or field == "code_multipath_f2":
                    
                for sys in sorted(self.dset.meta["obstypes"].keys()):
                    x_arrays = []
                    y_arrays = []
                    labels = []
                    solution = f"{field}_{sys}_{'_'.join(self.dset.meta['linear_combination'][field][sys])}"
                    figure_path = self.figure_dir / figure_name.replace("{solution}", solution).replace("{FIGURE_FORMAT}", self.figure_format)
                    figure_info.append(FigureInfo(
                                        figure_path, 
                                        field.replace("_f1", " ").replace("_f2", " "), 
                                        sys, 
                                        self.dset.meta['linear_combination'][field][sys]),
                    )
                        
                    for sat in sorted(self.dset.unique("satellite")):
                        if not sat.startswith(sys):
                            continue
                        idx = self.dset.filter(satellite= sat)
                        x_arrays.append(self.dset.time.gps.datetime[idx])
                        y_arrays.append(self.dset.lin[field][idx])
                        labels.append(sat)  

                    # Plot scatter plot
                    plt = MatPlotExt()
                    plt.plot(
                        x_arrays=x_arrays,
                        y_arrays=y_arrays,
                        xlabel="Time [GPS]",
                        ylabel="Code-multipath combination",
                        y_unit="m",
                        labels=labels,
                        figure_path=figure_path,
                        options={
                            "figsize": (7, 6),
                            "legend": True,
                            "legend_ncol": 6,
                            "legend_location": "bottom",
                            "plot_to": "file", 
                            "plot_type": "scatter",
                            "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
                            "xlim": "fit_to_data",
                        },
                    ) 
                    
            # Code-phase difference
            elif field == "code_phase_f1" or field == "code_phase_f2":
                    
                for sys in sorted(self.dset.meta["obstypes"].keys()):
                    x_arrays = []
                    y_arrays = []
                    labels = []
                    solution = f"{field}_{sys}_{'_'.join(self.dset.meta['linear_combination'][field][sys])}"
                    figure_path = self.figure_dir / figure_name.replace("{solution}", solution).replace("{FIGURE_FORMAT}", self.figure_format)
                    log.debug(f"Plot {figure_path}.")
                    figure_info.append(FigureInfo(
                                        figure_path, 
                                        field.replace("_f1", " ").replace("_f2", " "), 
                                        sys, 
                                        self.dset.meta['linear_combination'][field][sys]),
                    )
                        
                    for sat in sorted(self.dset.unique("satellite")):
                        if not sat.startswith(sys):
                            continue
                        idx = self.dset.filter(satellite= sat)
                        x_arrays.append(self.dset.time.gps.datetime[idx])
                        y_arrays.append(self.dset.lin[field][idx])
                        labels.append(sat)  

                    # Plot scatter plot
                    plt = MatPlotExt()
                    plt.plot(
                        x_arrays=x_arrays,
                        y_arrays=y_arrays,
                        xlabel="Time [GPS]",
                        ylabel="Code-phase difference",
                        y_unit="m",
                        labels=labels,
                        figure_path=figure_path,
                        options={
                            "figsize": (7, 6),
                            "legend": True,
                            "legend_ncol": 6,
                            "legend_location": "bottom",
                            "plot_to": "file", 
                            "plot_type": "scatter",
                            "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
                            "xlim": "fit_to_data",
                        },
                    ) 
                                        
            # Geometry_free, ionosphere_free, wide_lane and narrow_lane linear combination
            else:
                   
                for sys in sorted(self.dset.unique("system")):
                    x_arrays = []
                    y_arrays = []
                    labels = []
                    solution = f"{field}_{sys}_{'_'.join(self.dset.meta['linear_combination'][field][sys])}"
                    figure_path = self.figure_dir / figure_name.replace("{solution}", solution).replace("{FIGURE_FORMAT}", self.figure_format)
                    log.debug(f"Plot {figure_path}.")
                    figure_info.append(FigureInfo(
                                        figure_path, 
                                        "_".join(field.split('_')[0:2]), # Remove _code, _doppler, _phase or _snr  
                                        sys, 
                                        self.dset.meta['linear_combination'][field][sys]),
                    )
                        
                    for sat in sorted(self.dset.unique("satellite")):
                        if not sat.startswith(sys):
                            continue
                        idx = self.dset.filter(satellite= sat)
                        x_arrays.append(self.dset.time.gps.datetime[idx])
                        y_arrays.append(self.dset.lin[field][idx])
                        labels.append(sat)                        
                        
                    # Get ylabel
                    name1, name2, obscode = field.split("_")                                    
                    ylabel = f"{name1.capitalize()}-{name2} ({obscode})"
                    
                    # Plot scatter plot
                    plt = MatPlotExt()
                    plt.plot(
                        x_arrays=x_arrays,
                        y_arrays=y_arrays,
                        xlabel="Time [GPS]",
                        ylabel=ylabel,
                        y_unit="m",
                        labels=labels,
                        figure_path=figure_path,
                        options={
                            "figsize": (7, 6),
                            "legend": True,
                            "legend_ncol": 6,
                            "legend_location": "bottom",
                            "plot_to": "file", 
                            "plot_type": "scatter",
                            "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
                            "xlim": "fit_to_data",
                        },
                    )

        return figure_info


    def plot_number_of_satellites(
                        self, 
                        figure_name: str="plot_gnss_number_of_satellites_epoch.{FIGURE_FORMAT}",
    ) -> pathlib.PosixPath:
        """Plot number of satellites based for each GNSS
        
        Args:
            figure_name: File name of figure.
        
        Returns:
            Figure path.
        """
        
        # Generate x- and y-axis data per system
        x_arrays = []
        y_arrays = []
        labels = []
        figure_path = self.figure_dir / figure_name.replace("{FIGURE_FORMAT}", self.figure_format)
        log.debug(f"Plot {figure_path}.")
    
        for sys in sorted(self.dset.unique("system")):
            idx = self.dset.filter(system=sys)
            x_arrays.append(self.dset.time.gps.datetime[idx])
            y_arrays.append(gnss.get_number_of_satellites(
                                            self.dset.system[idx], 
                                            self.dset.satellite[idx], 
                                            self.dset.time.gps.datetime[idx])
            )
            labels.append(enums.gnss_id_to_name[sys].value)
    
        # Plot scatter plot
        plt = MatPlotExt()
        plt.plot(
            x_arrays=x_arrays,
            y_arrays=y_arrays,
            xlabel="Time [GPS]",
            ylabel="# satellites",
            y_unit="",
            labels=labels,
            figure_path=figure_path,
            options={
                "figsize": (7, 4), 
                "marker": ",",
                "legend": True,
                "legend_location": "bottom",
                "legend_ncol": len(self.dset.unique("system")),
                "plot_to": "file", 
                "plot_type": "plot",
                "xlim": "fit_to_data",
            },
        )
        
        return figure_path
    
    
    def plot_number_of_satellites_used(
            self, 
            figure_name: str="plot_number_of_satellites_used.{FIGURE_FORMAT}",
    ) -> pathlib.PosixPath:
        """Plot available number of satellites against used number
    
        Args:
            figure_name: File name of figure.
        
        Returns:
            Figure path.
        """
        figure_path = self.figure_dir / figure_name.replace("{FIGURE_FORMAT}", self.figure_format)
        log.debug(f"Plot {figure_path}.")
    
        if "num_satellite_used" not in self.dset.fields:
            self.dset.add_float(
                "num_satellite_used",
                val=gnss.get_number_of_satellites(
                                    self.dset.system, 
                                    self.dset.satellite, 
                                    self.dset.time.gps.datetime),
                write_level="operational",
            )
    
        plt = MatPlotExt()
        plt.plot(
            x_arrays=[self.dset.time.gps.datetime, self.dset.time.gps.datetime],
            y_arrays=[self.dset.num_satellite_available, self.dset.num_satellite_used],
            xlabel="Time [GPS]",
            ylabel="Number of satellites",
            y_unit="",
            labels=["Available", "Used"],
            figure_path=figure_path,
            options={
                "figsize": (7, 4), 
                "legend": True, 
                "marker": ",", 
                "plot_to": "file", 
                "plot_type": "plot",
                "xlim": "fit_to_data",
            },
        )
        
        return figure_path
    
        
    def plot_obstype_availability(
                            self, 
                            figure_name: str="plot_obstype_availability_{system}.{FIGURE_FORMAT}",
    ) -> List[pathlib.PosixPath]:
        """Generate GNSS observation type observation type availability based on RINEX observation file
 
        Args:
            figure_name: File name of figure.
            
        Returns:
            List with figure path for observation type availability depending on GNSS. File name ends with GNSS 
            identifier (e.g. 'E', 'G'), for example 'plot_obstype_availability_E.png'.
        """
        figure_paths = list()
    
        for sys in sorted(self.dset.unique("system")):

            x_arrays = []
            y_arrays = []
            
            idx_sys = self.dset.filter(system=sys)
            num_sat = len(set(self.dset.satellite[idx_sys]))
            
            figure_path = self.figure_dir / figure_name.replace("{system}", sys).replace("{FIGURE_FORMAT}", self.figure_format)
            figure_paths.append(figure_path)
            log.debug(f"Plot {figure_path}.")
            
            for sat in sorted(self.dset.unique("satellite"), reverse=True):
                if not sat.startswith(sys):
                    continue
                
                idx_sat = self.dset.filter(satellite=sat)
                keep_idx = np.full(self.dset.num_obs, False, dtype=bool)
                
                for obstype in self._sort_string_array(self.dset.meta["obstypes"][sys]):
                    keep_idx[idx_sat] = np.logical_not(np.isnan(self.dset.obs[obstype][idx_sat]))
                    #time_diff = np.diff(self.dset.time.gps.gps_seconds[keep_idx]) 
                    #time_diff = np.insert(time_diff, 0, float('nan')) 
                    
                    if np.any(keep_idx):
                        num_obs = len(self.dset.time[keep_idx])
                        x_arrays.append(self.dset.time.gps.datetime[keep_idx])
                        y_arrays.append(np.full(num_obs, f"{sat}_{obstype}"))

            plt = MatPlotExt()
            plt.plot(
                x_arrays=x_arrays,
                y_arrays=y_arrays,
                xlabel="Time [GPS]",
                ylabel="Satellite and observation type",
                figure_path=figure_path,
                y_unit="",
                options={
                    "colormap": "tab20",
                    "figsize": (1.0 * num_sat, 3.0 * num_sat),
                    "fontsize": 5,
                    "plot_to": "file",
                    "plot_type": "scatter",
                    #"title": "Satellite and observation type",
                    "xlim": "fit_to_data",
                },
            )

        return figure_paths
    

    def plot_satellite_availability(
                            self, 
                            figure_name: str="plot_satellite_availability.{FIGURE_FORMAT}",
    ) -> pathlib.PosixPath:
        """Generate GNSS satellite observation availability overview based on RINEX observation file
 
        Args:
            figure_name: File name of figure.
            
        Returns:
            Figure path.
        """
        
        # Generate x- and y-axis data per system
        x_arrays = []
        y_arrays = []
        labels = []
        figure_path = self.figure_dir / figure_name.replace("{FIGURE_FORMAT}", self.figure_format)
        log.debug(f"Plot {figure_path}.")
    
        time, satellite, system = self._sort_by_satellite()
    
        for sys in sorted(self.dset.unique("system"), reverse=True):
            idx = system == sys
            x_arrays.append(time[idx])
            y_arrays.append(satellite[idx])
            labels.append(enums.gnss_id_to_name[sys].value)
            
        # Plot scatter plot
        num_sat = len(self.dset.unique("satellite"))
        plt = MatPlotExt()
        plt.plot(
            x_arrays=x_arrays,
            y_arrays=y_arrays,
            xlabel="Time [GPS]",
            ylabel="Satellite",
            y_unit="",
            #labels=labels,
            figure_path=figure_path,
            options={
                "colormap": "tab20",
                "figsize": (0.1 * num_sat, 0.2 * num_sat),
                "fontsize": 10,
                "legend": True,
                "legend_location": "bottom",
                "legend_ncol": len(self.dset.unique("system")),
                "plot_to": "file",
                "plot_type": "scatter",
                #"title": "Satellite availability",
                "xlim": "fit_to_data",
            },
        )
        
        return figure_path
    
        
        
    def plot_skyplot(
                    self,
                    figure_name: str="plot_skyplot_{system}.{FIGURE_FORMAT}",
    ) -> List[pathlib.PosixPath]:
        """Plot skyplot for each GNSS
        
        Args:
            figure_name: File name of figure.
    
        Returns:
            List with figure path for skyplot depending on GNSS. File name ends with GNSS identifier (e.g. 'E', 'G'), 
            for example 'plot_skyplot_E.png'.
        """
        figure_paths = list()
    
        # Convert azimuth to range 0-360 degree
        try:
            azimuth = self.dset.site_pos.azimuth
        except exceptions.InitializationError:
            return figure_paths
        idx = azimuth < 0
        azimuth[idx] = 2 * np.pi + azimuth[idx]
    
        # Convert zenith distance from radian to degree
        zenith_distance = np.rad2deg(self.dset.site_pos.zenith_distance)
    
        # Generate x- and y-axis data per system
        for sys in sorted(self.dset.unique("system")):
            x_arrays = []
            y_arrays = []
            labels = []
           
            figure_path = self.figure_dir / figure_name.replace("{system}", sys).replace("{FIGURE_FORMAT}", self.figure_format)
            figure_paths.append(figure_path)
            
            for sat in sorted(self.dset.unique("satellite")):
                if not sat.startswith(sys):
                    continue
                idx = self.dset.filter(satellite= sat)
                x_arrays.append(azimuth[idx])
                y_arrays.append(zenith_distance[idx])
                labels.append(sat)
        
            # Plot with polar projection
            # TODO: y-axis labels are overwritten after second array plot. Why? What to do?
            plt = MatPlotExt()
            plt.plot(
                x_arrays=x_arrays,
                y_arrays=y_arrays,
                xlabel="",
                ylabel="",
                y_unit="",
                labels=labels,
                figure_path=figure_path,
                options={
                    "colormap": "tab20",
                    "figsize": (7, 7.5),
                    "legend": True,
                    "legend_ncol": 6,
                    "legend_location": "bottom",
                    "plot_to": "file",
                    "plot_type": "scatter",
                    "projection": "polar",
                    "title": f"Skyplot for {enums.gnss_id_to_name[sys]}\n Azimuth [deg] / Elevation[deg]",
                    "xlim": [0, 2 * np.pi],
                    "ylim": [0, 90],
                    "yticks": (range(0, 90, 30)),  # sets 3 concentric circles
                    "yticklabels": (map(str, range(90, 0, -30))),  # reverse labels from zenith distance to elevation
                },
            )
        
        return figure_paths
    
    
    def plot_satellite_elevation(
                    self,
                    figure_name: str="plot_satellite_elevation_{system}.{FIGURE_FORMAT}",
    ) -> List[pathlib.PosixPath]:
        """Plot satellite elevation for each GNSS
    
        Args:
            figure_name: File name of figure.
            
        Returns:
            List with figure path for skyplot depending on GNSS. File name ends with GNSS identifier (e.g. 'E', 'G'), 
            for example 'plot_skyplot_E.png'.
        """
        figure_paths = list()
 
        # Convert elevation from radian to degree
        try:
            elevation = np.rad2deg(self.dset.site_pos.elevation)
        except exceptions.InitializationError:
            return figure_paths
            
    
        # Limit x-axis range to rundate
        day_start, day_end = get_day_limits(self.dset)
         
        # Generate x- and y-axis data per system
        for sys in sorted(self.dset.unique("system")):
            x_arrays = []
            y_arrays = []
            labels = []
            
            figure_path = self.figure_dir / figure_name.replace("{system}", sys).replace("{FIGURE_FORMAT}", self.figure_format)
            figure_paths.append(figure_path)
            log.debug(f"Plot {figure_path}.")
    
            for sat in sorted(self.dset.unique("satellite")):
                if not sat.startswith(sys):
                    continue
                idx = self.dset.filter(satellite=sat)
                x_arrays.append(self.dset.time.gps.datetime[idx])
                y_arrays.append(elevation[idx])
                labels.append(sat)
        
            # Plot with scatter plot
            plt = MatPlotExt()
            plt.plot(
                x_arrays=x_arrays,
                y_arrays=y_arrays,
                xlabel="Time [GPS]",
                ylabel="Elevation [deg]",
                y_unit="",
                labels=labels,
                figure_path=figure_path,
                options={
                    "colormap": "tab20",
                    "figsize": (7, 8),
                    "legend": True,
                    "legend_ncol": 6,
                    "legend_location": "bottom",
                    "plot_to": "file",
                    "plot_type": "scatter",
                    "title": f"Satellite elevation for {enums.gnss_id_to_name[sys]}",
                    "xlim": [day_start, day_end],
                },
            )
            
        return figure_paths



    def plot_satellite_overview(
                    self,
                    figure_name: str="plot_satellite_overview.{FIGURE_FORMAT}",
    ) -> Union[pathlib.PosixPath, None]:
        """Plot satellite observation overview
    
        Args:
           dset:        A dataset containing the data.
           figure_dir:  Figure directory
           
        Returns:
           Figure path or None if necessary datasets could not be read
        """
        figure_path = self.figure_dir / figure_name.replace("{FIGURE_FORMAT}", self.figure_format)
        log.debug(f"Plot {figure_path}.")
    
        # Limit x-axis range to rundate
        day_start, day_end = get_day_limits(self.dset)
    
        # Get time and satellite data from read and orbit stage
        file_vars = {**self.dset.vars, **self.dset.analysis}
        file_vars["stage"] = "read"
        file_path = config.files.path("dataset", file_vars=file_vars)
        if file_path.exists():
            systems = list(self.dset.unique("system")) # self.dset.meta["obstypes"].keys()
            time_read, satellite_read, _ = self._sort_by_satellite(
                self._get_dataset(stage="read", systems=systems)
            )
            time_orbit, satellite_orbit, _ = self._sort_by_satellite(
                self._get_dataset(stage="orbit", systems=systems)
            )
            time_edit, satellite_edit, _ = self._sort_by_satellite(
                self._get_dataset(stage="edit", systems=systems)
            )
            
        else:
            # NOTE: This is the case for concatencated Datasets, where "read" and "edit" stage data are not available.
            log.warn(f"Read dataset does not exists: {file_path}. Plot {figure_path} can not be plotted.")
            return None
    
        # Generate plot
        plt = MatPlotExt()
        plt.plot(
            x_arrays=[time_read.tolist(), time_orbit.tolist(), time_edit.tolist()],
            y_arrays=[satellite_read.tolist(), satellite_orbit.tolist(), satellite_edit.tolist()],
            xlabel="Time [GPS]",
            ylabel="Satellite",
            y_unit="",
            # labels = ["Rejected in orbit stage", "Rejected in edit stage", "Kept observations"],
            colors=["red", "orange", "green"],
            figure_path=figure_path,
            options={
                "colormap": "tab20",
                "figsize": (7, 6),
                "marker": "|",
                "plot_to": "file",
                "plot_type": "scatter",
                "title": "Overview over satellites",
                "xlim": [day_start, day_end],
            },
        )
        
        return figure_path
    
        
    def plot_tgd_comparison(
                            self, 
                            figure_name: str="plot_{solution}.{FIGURE_FORMAT}",

    ) -> List[pathlib.PosixPath]:
        """Plot total/broadcast group delay (TGD/BGD) comparison results
        
        Args:
            figure_name: File name of figure.
        
        Returns:
            List with figure path for TGD/BGD comparison plots. File name ends with GNSS identifier
            (e.g. 'E', 'G') and observation type, for example 'plot_G_tgd.png'.
        """          
        figure_paths = list()
       
        plot_def = {
            "C": {
                "tgd_b1_b3": PlotField("TGD(B1,B2), DCB(C2I,C6I)", "ns", (7, 4)),
                "tgd_b2_b3": PlotField("TGD(B1,B3), DCB(C7I,C6I)", "ns", (7, 4)),
            },
            "E": {
                "bgd_e1_e5a":  PlotField("BGD(E1,E5a), DCB(C1C,C5Q)", "ns", (7, 4)),
                "bgd_e1_e5b":  PlotField("BGD(E1,E5b), DCB(C1C,C7Q)", "ns", (7, 4)),
            },
            "G": {
                "tgd": PlotField("TGD(L1,L2), DCB(C1W,C2W)", "ns", (7, 4)),
            },
            "J": {
                "tgd": PlotField("TGD(L1,L2), DCB(C1X,C2X)", "ns", (7, 4)),
            },
           
        }  
            
        tgd_def = ["bgd_e1_e5a", "bgd_e1_e5b", "tgd", "tgd_b1_b3", "tgd_b2_b3"]

        #
        # Plot average of TGD and DCB by satellite
        #          

        # Generate x- and y-axis data per system
        for sys in sorted(self.dset.unique("system")):
             
            if sys not in plot_def:
                log.warn(f"TGD/BGD comparison against post-processed DCBs is not defined for GNSS '{sys}'.")
                continue
          
            for field in plot_def[sys].keys():

                if not field in self.dset.fields or not f"{field}_dcb" in self.dset.fields:
                    continue
            
                satellites = []
                tgd_sat_mean = []
                dcb_sat_mean = []
                tgd_dcb_sat_mean = []
                tgd_dcb_sat_std = []
                tgd_dcb_sat_rms = []
                tgd_dcb_sat_percentile = []
                        
                for sat in sorted(self.dset.unique("satellite")):
                    if not sat.startswith(sys):
                        continue
                    idx = self.dset.filter(satellite=sat)
                    satellites.append(sat)
                    tgd_sat = self.dset[f"{field}_mean"][idx] * Unit.second2nanosecond
                    dcb_sat = self.dset[f"{field}_dcb_mean"][idx] * Unit.second2nanosecond
                    tgd_dcb_sat = tgd_sat - dcb_sat
                    tgd_sat_mean.append(np.mean(tgd_sat))
                    dcb_sat_mean.append(np.mean(dcb_sat))
                    tgd_dcb_sat_mean.append(np.mean(tgd_dcb_sat))
                    tgd_dcb_sat_std.append(np.std(tgd_dcb_sat))
                    tgd_dcb_sat_rms.append(self._rms(tgd_dcb_sat))
                    tgd_dcb_sat_percentile.append(np.nanpercentile((np.absolute(tgd_dcb_sat)), 95))
            
                # Plot average of TGD and DCB by satellite in one plot
                figure_path = self.figure_dir / figure_name.replace("{FIGURE_FORMAT}", self.figure_format).replace("{solution}", f"{sys}_{field}_dcb")
                figure_paths.append(figure_path)
                log.debug(f"Plot {figure_path}.")
                         
                plt = MatPlotExt()
                plt.plot(
                    x_arrays= [satellites, satellites],
                    y_arrays=[tgd_sat_mean, dcb_sat_mean],
                    xlabel="Satellite",
                    ylabel=f"Mean {plot_def[sys][field].ylabel}",
                    y_unit="ns",
                    labels=plot_def[sys][field].ylabel.split(", "),
                    colors=["blue", "red"],
                    figure_path=figure_path,
                    options={
                        "figsize": (7, 4),
                        "legend": True,
                        "legend_ncol": 2,
                        "legend_location": "bottom",
                        "marker": ".",
                        "markersize": 15,
                        "plot_to": "file",
                        "plot_type": "scatter",
                        "xlabelrotation": 45,
                    },
                )
                               
                # Plot difference average of TGD and DCB by satellite 
                figure_path = self.figure_dir / figure_name.replace("{FIGURE_FORMAT}", self.figure_format).replace("{solution}", f"{sys}_{field}_dcb_diff")
                figure_paths.append(figure_path)
                log.debug(f"Plot {figure_path}.")
                plt = MatPlotExt()
                plt.plot(
                    x_arrays= [satellites],
                    y_arrays=[tgd_dcb_sat_mean],
                    yerr_arrays=[tgd_dcb_sat_std],
                    xlabel="Satellite",
                    ylabel=f"Mean {plot_def[sys][field].ylabel.replace(', ', ' - ')}",
                    y_unit="ns",
                    colors=["black"],
                    figure_path=figure_path,
                    options={
                        "ecolor": "grey",
                        "errorbar": True,
                        "figsize": plot_def[sys][field].figsize,
                        "marker": ".",
                        "markersize": 30,
                        "plot_to": "file",
                        "plot_type": "scatter",
                        "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
                        "xlabelrotation": 45,
                    },
                )
                
                # RMS
                figure_path = self.figure_dir / figure_name.replace("{FIGURE_FORMAT}", self.figure_format).replace("{solution}", f"{sys}_{field}_dcb_diff_rms")
                figure_paths.append(figure_path)
                log.debug(f"Plot {figure_path}.")
                plt = MatPlotExt()
                plt.plot(
                    x_arrays= [satellites],
                    y_arrays=[tgd_dcb_sat_rms],
                    xlabel="Satellite",
                    ylabel=f"RMS {plot_def[sys][field].ylabel.replace(', ', ' - ')}",
                    y_unit="ns",
                    colors=["dodgerblue"],
                    figure_path=figure_path,
                    options={
                        "figsize": plot_def[sys][field].figsize,
                        "marker": ".",
                        "markersize": 30,
                        "plot_to": "file",
                        "plot_type": "bar",
                        "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
                        "xlabelrotation": 45,
                    },
                )
                
                # Percentile
                figure_path = self.figure_dir / figure_name.replace("{FIGURE_FORMAT}", self.figure_format).replace("{solution}", f"{sys}_{field}_dcb_diff_percentile")
                figure_paths.append(figure_path)
                log.debug(f"Plot {figure_path}.")
                plt = MatPlotExt()
                plt.plot(
                    x_arrays= [satellites],
                    y_arrays=[tgd_dcb_sat_percentile],
                    xlabel="Satellite",
                    ylabel=f"Percentile {plot_def[sys][field].ylabel.replace(', ', ' - ')}",
                    y_unit="ns",
                    colors=["dodgerblue"],
                    figure_path=figure_path,
                    options={
                        "figsize": plot_def[sys][field].figsize,
                        "marker": ".",
                        "markersize": 30,
                        "plot_to": "file",
                        "plot_type": "bar",
                        "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
                        "xlabelrotation": 45,
                    },
                )
                

        #
        # Plot TGD, TGD-DCB and TGD_mean - DCB_mean by time
        #         
        for field in tgd_def: 
            
            if field in self.dset.fields:
                figure_paths = figure_paths + self.plot_field(field)
            
            #if f"{field}_diff" in self.dset.fields:
            #    figure_paths = figure_paths + self.plot_field(f"{field}_diff")
                
            if f"{field}_diff_mean" in self.dset.fields:
                figure_paths = figure_paths + self.plot_field(f"{field}_diff_mean")
        
        return figure_paths


    #        
    # ADDITIONAL PUBLIC FUNCTIONS
    # 
    def get_first_galileo_signal(self):
        """Get first Galileo signal given in the navigation message based on ordered list
       
        Returns:
           Tuple with chosen Galileo signal and navigation message type
        """        
        signals = self._select_gnss_signal(system="E")
        signal = next(iter(sorted(signals)))
        nav_type = signals[signal]
    
        return signal, nav_type
    
    
    #        
    # AUXILIARY FUNCTIONS
    #
    def _convert_to_unit(self, data: np.ndarray, from_unit: str, to_unit: str) -> np.ndarray:
        """Convert field data to specified unit
        
        Args:
            data:       Array with data
            from_unit:  Data unit
            to_unit:    Unit to convert 
            
        Returns:
            Array with data in specified unit
        """
        if (from_unit is None) or (from_unit == to_unit):
            return data

        #log.debug(f"Convert field {field} from unit {from_unit} to {to_unit}.")
        return data * Unit(from_unit).to(to_unit).m 
    
    
    def _get_dataset(
            self, 
            stage: str, 
            systems: Union[List[str], None] = None,
    ) -> "Dataset":
        """Get dataset for given stage
    
        Args:
           systems:     List with GNSS identifiers (e.g. E, G, ...)
    
        Returns:
           Dataset for given stage or error exit status if dataset could not be read
        """
    
        # Get Dataset
        # TODO: "label" should have a default value.
        file_vars = {**self.dset.vars, **self.dset.analysis}
        file_vars["stage"] = stage
        try:
            dset_out = dataset.Dataset.read(**file_vars)
        except OSError:
            log.warn("Could not read dataset {config.files.path('dataset', file_vars=file_vars)}.")
            return enums.ExitStatus.error
    
        # Reject not defined GNSS observations
        if systems:
            systems = [systems] if isinstance(systems, str) else systems
            keep_idx = np.zeros(dset_out.num_obs, dtype=bool)
            for sys in systems:
                idx = dset_out.filter(system=sys)
                keep_idx[idx] = True
            dset_out.subset(keep_idx)
    
        return dset_out
    
    

    
        
    def _get_gnss_signal_in_space_status_data(
                self,
                status: int, 
                signal: Union[str, None]=None,
                system: Union[str, None]=None,
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
           status:        Signal in space status.   
           signal:        Galileo signal for which SIS status should be determined. The signal can be e1, e5a or e5b.
           system:        Get only signal in space status for selected system identifier.
    
        Returns:
            Tuple with time and satellite data for a given SIS status and ordered by satellite
        """
        # TODO: order of satellites is not correct. Maybe to save the data in a dataframe could help?
    
        # Generate x- and y-axis data
        time = []
        satellite = []
    
        for sat in sorted(self.dset.unique("satellite"), reverse=True):
            
            if system: # Skip GNSS, which are not selected
                if not sat.startswith(system):
                    continue
                
            idx = self.dset.filter(satellite=sat)
    
            if sat.startswith("E"):
                idx_status = getattr(self.dset, "sis_status_" + signal)[idx] == status
            else:
                if status == 0:  # healthy status
                    idx_status = self.dset.sv_health[idx] == 0
                elif status == 3:  # unhealthy status
                    idx_status = self.dset.sv_health[idx] > 0
                else:
                    continue
    
            time.extend(self.dset.time.gps.datetime[idx][idx_status])
            satellite.extend(self.dset.satellite[idx][idx_status])
    
        return time, satellite
    
        
    def _rms(self, x: np.ndarray) -> np.ndarray:
        """Determine root mean square (RMS)
        
        Args:
            x: Array with data
            
        Returns:
            RMS for given data
        """
        return np.sqrt(np.nanmean(np.square(x)))
        
        
    def _select_gnss_signal(self, system: Union[str, None]=None) -> Tuple[List[str], str]:
        """Select GNSS signal depending on given data in Dataset
    
        Args:
            system:  System identifier
            
        Returns:
            Selected GNSS signal
        """
        nav_type_def = {
                "FNAV_E5a": {"e5a": "F/NAV"},
                "INAV_E1": {"e1": "I/NAV"},
                "INAV_E5b": {"e5b": "I/NAV"},
                "INAV_E1E5b": {"e1": "I/NAV", "e5b": "I/NAV"},
                "LNAV": {None: "L/NAV"},
                "D1/D2": {None: "D1/D2"},
                "NAV": {None: "NAV"},
        }
        
        signals = dict()
        
        if system:
            idx = self.dset.filter(system=system)
            nav_types = set(self.dset.nav_type[idx])
        else:
            nav_types = self.dset.unique("nav_type")
            
        for nav_type in nav_types:
            try:
                signals.update(nav_type_def[nav_type])
            except KeyError:
                log.error(f"GNSS navigation message {nav_type} is not defined.")
                
        return signals
    
      
    def _sort_by_satellite(self, dset: Union["Dataset", None]=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Sort time and satellite fields of dataset by satellite order
    
        Returns:
            Tuple with ordered time, satellite and system array
        """
        dset = self.dset if dset is None else dset
        time = []
        satellite = []
        system = []
        for sat in sorted(self.dset.unique("satellite"), reverse=True):
            idx = self.dset.filter(satellite=sat)
            time.extend(self.dset.time.gps.datetime[idx])
            satellite.extend(self.dset.satellite[idx])
            system.extend(self.dset.system[idx])
    
        return np.array(time), np.array(satellite), np.array(system)
    
    
    def _sort_string_array(self, array: List[str]) -> List[str]:
        """Sort string array based on last two characters
    
        Args:
            array: String array
            
        Returns:
            Sorted string array
        """
        array.sort(key=lambda x: x[-2:3])
        return array
