"""Class for GNSS plots

Description:
------------
TODO

"""
# Standard liberay imports
from datetime import datetime
import pathlib
from typing import List, Tuple

# External liberary imports
import numpy as np

# Midgard imports
from midgard.collections import enums
from midgard.gnss import gnss
from midgard.plot.matplotlib_extension import plot


FIGURE_FORMAT = "png"

class GnssPlot:
    """Class for GNSS plots
    """

    def __init__(
        self, 
        dset: "Dataset",
        figure_dir: pathlib.PosixPath,
    ) -> None:
        """Set up a new GNSS plot object

        Args:
            dset:          A dataset containing the data.
            figure_dir:    Figure directory.
        """
        self.dset = dset
        self.figure_dir = figure_dir


    def plot_number_of_satellites(
                        self, 
                        figure_name: str=f"plot_gnss_number_of_satellites_epoch.{FIGURE_FORMAT}",
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
        figure_path = self.figure_dir / figure_name
    
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
        plot(
            x_arrays=x_arrays,
            y_arrays=y_arrays,
            xlabel="Time [GPS]",
            ylabel="# satellites",
            y_unit="",
            labels=labels,
            figure_path=figure_path,
            opt_args={
                "figsize": (7, 4), 
                "marker": ",",
                "legend": True,
                "legend_location": "bottom",
                "legend_ncol": len(self.dset.unique("system")),
                "plot_to": "file", 
                "plot_type": "plot"
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
            
            figure_path = self.figure_dir / figure_name.replace("{system}", sys).replace("{FIGURE_FORMAT}", FIGURE_FORMAT)
            figure_paths.append(figure_path)
            
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

            plot(
                x_arrays=x_arrays,
                y_arrays=y_arrays,
                xlabel="Time [GPS]",
                ylabel="Satellite and observation type",
                figure_path=figure_path,
                y_unit="",
                opt_args={
                    "colormap": "tab20",
                    "figsize": (1.0 * num_sat, 3.0 * num_sat),
                    "fontsize": 5,
                    "plot_to": "file",
                    "plot_type": "scatter",
                    #"title": "Satellite and observation type",
                },
            )

        return figure_paths
    

    def plot_satellite_availability(
                            self, 
                            figure_name: str=f"plot_satellite_availability.{FIGURE_FORMAT}",
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
        figure_path = self.figure_dir / figure_name
    
        time, satellite, system = self._sort_by_satellite()
    
        for sys in sorted(self.dset.unique("system"), reverse=True):
            idx = system == sys
            x_arrays.append(time[idx])
            y_arrays.append(satellite[idx])
            labels.append(enums.gnss_id_to_name[sys].value)
            
        # Plot scatter plot
        num_sat = len(self.dset.unique("satellite"))
        plot(
            x_arrays=x_arrays,
            y_arrays=y_arrays,
            xlabel="Time [GPS]",
            ylabel="Satellite",
            y_unit="",
            #labels=labels,
            figure_path=figure_path,
            opt_args={
                "colormap": "tab20",
                "figsize": (0.1 * num_sat, 0.2 * num_sat),
                "fontsize": 10,
                "legend": True,
                "legend_location": "bottom",
                "legend_ncol": len(self.dset.unique("system")),
                "plot_to": "file",
                "plot_type": "scatter",
                #"title": "Satellite availability",
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
        azimuth = self.dset.site_pos.azimuth
        idx = azimuth < 0
        azimuth[idx] = 2 * np.pi + azimuth[idx]
    
        # Convert zenith distance from radian to degree
        zenith_distance = np.rad2deg(self.dset.site_pos.zenith_distance)
    
        # Generate x- and y-axis data per system
        for sys in sorted(self.dset.unique("system")):
            x_arrays = []
            y_arrays = []
            labels = []
           
            figure_path = self.figure_dir / figure_name.replace("{system}", sys).replace("{FIGURE_FORMAT}", FIGURE_FORMAT)
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
            plot(
                x_arrays=x_arrays,
                y_arrays=y_arrays,
                xlabel="",
                ylabel="",
                y_unit="",
                labels=labels,
                figure_path=figure_path,
                opt_args={
                    "colormap": "hsv",
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
        elevation = np.rad2deg(self.dset.site_pos.elevation)
    
        # Limit x-axis range to rundate
        day_start, day_end = self._get_day_limits()
         
        # Generate x- and y-axis data per system
        for sys in sorted(self.dset.unique("system")):
            x_arrays = []
            y_arrays = []
            labels = []
            
            figure_path = self.figure_dir / figure_name.replace("{system}", sys).replace("{FIGURE_FORMAT}", FIGURE_FORMAT)
            figure_paths.append(figure_path)
    
            for sat in sorted(self.dset.unique("satellite")):
                if not sat.startswith(sys):
                    continue
                idx = self.dset.filter(satellite=sat)
                x_arrays.append(self.dset.time.gps.datetime[idx])
                y_arrays.append(elevation[idx])
                labels.append(sat)
        
            # Plot with scatter plot
            plot(
                x_arrays=x_arrays,
                y_arrays=y_arrays,
                xlabel="Time [GPS]",
                ylabel="Elevation [deg]",
                y_unit="",
                labels=labels,
                figure_path=figure_path,
                opt_args={
                    "colormap": "hsv",
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


    #        
    # AUXILIARY FUNCTIONS
    # 
    def _get_day_limits(self) -> Tuple[datetime, datetime]:
        """Get start and end time for given run date
    
            Returns:
                Start and end date. 
            """
        day_start = min(self.dset.time.datetime)
        day_end = max(self.dset.time.datetime)
    
        return day_start, day_end

      
    def _sort_by_satellite(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Sort time and satellite fields of dataset by satellite order
    
        Returns:
            Tuple with ordered time, satellite and system array
        """
        time = []
        satellite = []
        system = []
        for sat in sorted(self.dset.unique("satellite"), reverse=True):
            idx = self.dset.filter(satellite=sat)
            time.extend(self.dset.time.gps.datetime[idx])
            satellite.extend(self.dset.satellite[idx])
            system.extend(self.dset.system[idx])
    
        return np.array([time]), np.array([satellite]), np.array([system])
    
    
    def _sort_string_array(self, array: List[str]) -> List[str]:
        """Sort string array based on last two characters
    
        Args:
            array: String array
            
        Returns:
            Sorted string array
        """
        array.sort(key=lambda x: x[-2:3])
        return array