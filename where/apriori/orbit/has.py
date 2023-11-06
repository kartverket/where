"""Handling of Galileo high accuracy service (HAS) orbits

Description:
------------
The module includes a class for handling apriori Galileo high accuracy service (HAS) orbits.

Example:
    from where import apriori

    # Get HasOrbit object
    has = apriori.get('orbit', rundate=rundate, station=station, system=system, apriori_orbit='has')

    # Write calculated Dataset to file
    has.dset.write()

"""

# Standard library imports
from datetime import timedelta
from typing import List, Union

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math import nputil

# Where imports
from where import apriori, cleaners, parsers
from where.apriori import orbit
from where.lib import config, log, util


@plugins.register
class HasOrbit(orbit.AprioriOrbit):
    """A class for representing Galileo high accuracy service (HAS) orbits.

    HAS corrections based on JRC HAS decoder can be read and edited. In addition HAS satellite position and satellite clock
    correction can be calculated for given observation epochs based on broadcast ephemeris.

    Attributes:
        dset (Dataset):         Dataset object, which includes satellite position and velocity, satellite clock
                                correction and relativistic clock correction for each observation epoch
        dset_edit (Dataset):    Dataset object, which includes edited broadcast ephemeris navigation messages
        dset_raw (Dataset):     Dataset object, which includes broadcast ephemeris navigation messages read from RINEX
                                file
        file_key (str):         Key to the broadcast orbit file defined in files.conf file.
        name (str):             Apriori orbit name
        system (tuple):         GNSS system identifiers

    Methods:
        satellite_clock_correction():  Determine satellite clock correction
        _calculate():           Calculate HAS satellite orbit and satellite clock correction for given observation
                                epochs
        _edit():                Edit JRC HAS decoder data and save it in a Dataset
        _read():                Read JRC HAS decoder data and save it in a Dataset
    """

    name = "has"

    def __init__(self, rundate, file_key: str, day_offset: int=1) -> None:
        """Set up a new HasOrbit object, does not parse any data

        TODO: Remove dependency on rundate, use time to read correct files. (What to do with dataset?)

        Args:
            rundate (date):     Date of model run.
            file_key (str):     File key of Galileo HAS file.
            day_offset (int):   Day offset used to calculate the number of days to read.
        """
        super().__init__(rundate=rundate)
        self.file_key = file_key
        self.day_offset = day_offset
        
        
    def _read(self, dset_raw):
        """Read JRC HAS decoder messages depending on given `file_key`

        Note that beside the given day also the HAS messages files from the day/days before and the day/days after can
        be read depending on given `day_offset`.

        Args:
            dset_raw:     Dataset representing raw data from HAS messages
        """
        date_to_read = dset_raw.analysis["rundate"] - timedelta(days=self.day_offset)
        file_paths = list()

        # Loop over days to read
        while date_to_read <= dset_raw.analysis["rundate"] + timedelta(days=self.day_offset):

            date = date_to_read.strftime("%Y-%m-%d")
            meta = dict()
            
            file_path = config.files.path(self.file_key, file_vars={**dset_raw.vars, **dset_raw.analysis})
            if not file_path.exists():                       
                log.fatal(f"File does not exists: {file_path}")
            
            if util.is_file_empty(file_path):
                log.fatal(f"File is empty: {file_path}")
                
            p = parsers.parse("gnss_has_decoder", file_path=file_path)
            if not p.data_available:
                log.fatal(f"No observations in file {file_path}.")
                
            file_paths.append(str(file_path))
            dset_temp = p.as_dataset()          

            # Extend Dataset dset_raw with temporary Dataset
            if dset_raw.num_obs:
                dset_raw.extend(dset_temp, meta_key=date)

            else:
                # Add date key to meta TODO: Better solution?
                meta.setdefault(date, dict()).update(dset_temp.meta)
                for k in dict(dset_temp["meta"]).keys():
                    del dset_temp.meta[k]
                dset_temp.meta.update(meta)

                dset_raw.update_from(p.as_dataset())

            date_to_read += timedelta(days=1)

        dset_raw.meta.add("file_path", file_paths, section="parser")
        
        
    def _edit(self, dset_edit):
        """Edit HAS message file data 
        
        Remove data from HAS messages, which are declared as not available or which should not be used. This is the
        case for:
            HAS orbit data    - Fields `delta_cross_track`, `delta_in_track`, delta_radial` are set to Not a number (NaN)
            HAS clock data    - `status` field are set to 1 or 2
            HAS bias data     - `av_flag` field is set to 0
        
        TODO: Is that necessary? First the navigation messages are sorted after the satellite and time of transmission.
              Removing of duplicated HAS message records.

        Args:
            dset_edit (Dataset):     Dataset representing edited data from RINEX navigation file
        """
        
        # Edit HAS orbit message data
        if "delta_cross_track" in self.dset_raw.fields:
            remove_idx = np.isnan(self.dset_raw.delta_cross_track)            
            log.info(f"Removing {sum(np.logical_not(~remove_idx)):5d} not available HAS orbit message data (declared as NaN).")

        # Edit HAS clock message data        
        if "status" in self.dset_raw.fields:
            remove_idx = self.dset_raw.status > 0
            log.info(f"Removing {sum(np.logical_not(~remove_idx)):5d} not available HAS clock message data (status=1) " 
                     f"and data which shall not be used (status=2).")
            
        # Edit HAS code and phase bias message data        
        if "av_flag" in self.dset_raw.fields:
            remove_idx = self.dset_raw.av_flag == 0
            log.info(f"Removing {sum(np.logical_not(~remove_idx)):5d} not available HAS code/phase bias message data (av_flag=0)." )
            
        # Copy raw dataset
        dset_edit.update_from(self.dset_raw)

        log.info(f"Keeping {sum(~remove_idx)} of {dset_edit.num_obs} HAS messages.")
        dset_edit.subset(~remove_idx)
        
        if dset_edit.num_obs == 0:
            log.fatal("No observations are available.")

           
    def add_orbit_to_dataset(self, dset: "Dataset") -> None:
        """Add HAS orbit correction to dataset 
        
        Args: 
            dset:   A Dataset containing HAS message data.
        """
        # Clean orbits by removing unavailable satellites, unhealthy satellites and checking validity length of
        # navigation records
        cleaners.apply_remover("gnss_clean_orbit_has", dset, orbit=self)

        # Get correct HAS message for given observations times by determining the indices to HAS message dataset
        dset_has_idx = self._get_has_message_idx(dset)
        
        # Get reference IOD of HAS orbit correction, which is needed to select correct broadcast navigation message
        dset.add_float(
            "has_gnssiod_orb", 
            val=np.take(self.dset_edit.gnssiod, dset_has_idx), 
            unit="meters", 
            write_level="operational",
        )
        
        # Add additional information       
        dset.add_float("has_iod", val=np.take(self.dset_edit.iod, dset_has_idx), unit="meters", write_level="analysis")
        dset.add_float(
            "has_validity_orb", 
            val=np.take(self.dset_edit.validity, dset_has_idx),
            unit="",
            write_level="analysis",
        )
        dset.add_float(
            "has_delta_radial", 
            val=np.take(self.dset_edit.delta_radial, dset_has_idx),
            unit="meter",
            write_level="analysis",
        )
        dset.add_float(
            "has_delta_in_track", 
            val=np.take(self.dset_edit.delta_in_track, dset_has_idx),
            unit="meter",
            write_level="analysis",
        )
        dset.add_float(
            "has_delta_cross_track", 
            val=np.take(self.dset_edit.delta_cross_track, dset_has_idx),
            unit="meter",
            write_level="analysis",
        )       
        dset.add_time(
            "has_reception_time_of_message_orb", 
            val=np.take(self.dset_edit.time.gps.datetime, dset_has_idx),
            scale="gps",
            fmt="datetime",
            write_level="analysis",
        )
        dset.add_time(
            "has_time_of_message_orb", 
            val=np.take(self.dset_edit.tom.gps.datetime, dset_has_idx),
            scale="gps",
            fmt="datetime",
            write_level="analysis",
        )
                
        # Determine orbits for corresponding reference IOD of HAS messages
        # Get broadcast orbit and clocks
        brdc = apriori.get(
            "orbit", 
            rundate=dset.analysis["rundate"], 
            system=tuple(dset.unique("system")), 
            station=dset.vars["station"],
            day_offset=1,  # This guarantees, that corresponding broadcast messages can be found especially at the 
                           # start of the day.
            apriori_orbit="broadcast",
        )
        brdc.calculate_orbit(dset, time="time")        
        dset.add_posvel(
            "sat_posvel", 
            time=dset.time, 
            system="trs", 
            val=brdc.dset.sat_posvel.trs, 
            write_level="operational",
        )
        dset.add_float(
            "delay.gnss_satellite_clock", 
            val=-brdc.dset.gnss_satellite_clock, 
            unit="meter",
            write_level="analysis",
        )
        dset.add_float(
            "delay.gnss_relativistic_clock", 
            val=-brdc.dset.gnss_relativistic_clock, 
            unit="meter",
            write_level="analysis",
        )
        
        dset.add_float("used_iode", val=brdc.dset.used_iode, write_level="analysis")
        dset.add_time("used_toe", val=brdc.dset.used_toe, write_level="analysis")
        dset.add_time("used_transmission_time", val=brdc.dset.used_transmission_time, write_level="analysis")
        
        # Determine HAS orbit correction 
        dhas_ric = np.stack([
                        dset.has_delta_radial, 
                        dset.has_delta_in_track, 
                        dset.has_delta_cross_track,   
        ], axis=1)
        
        dhas_trs = np.squeeze(self._ric2trs(dset.sat_posvel) @ dhas_ric[:, :, None])

        dset.add_posvel_delta(
            name="has_orbit_correction",
            val=np.hstack([dhas_trs, np.zeros((dset.num_obs, 3))]),
            system="trs",
            ref_pos=dset.sat_posvel,
            write_level="operational",
        )
        
                     
    def add_clock_to_dataset(self, dset: "Dataset") -> None:
        """Add HAS clock correction to dataset
        
        Args: 
            dset:   A Dataset containing model data.
        """
        # Clean orbits by removing unavailable satellites, unhealthy satellites and checking validity length of
        # navigation records
        cleaners.apply_remover("gnss_clean_orbit_has", dset, orbit=self)

        # Get correct HAS message for given observations times by determining the indices to HAS message dataset
        dset_has_idx = self._get_has_message_idx(dset)
        
        # Add additional information
        dset.add_float(
            "has_gnssiod_clk", 
            val=np.take(self.dset_edit.gnssiod, dset_has_idx), 
            unit="meters", 
            write_level="operational",
        )
        
        dset.add_float(
            "has_clock_correction", 
            val=self.dset_edit.delta_clock_c0[dset_has_idx] * self.dset_edit.multiplier[dset_has_idx], 
            unit="meters",
            write_level="analysis",
        )
        dset.add_float(
            "has_validity_clk", 
            val=np.take(self.dset_edit.validity, dset_has_idx),
            unit="",
            write_level="analysis",
        )
        dset.add_time(
            "has_reception_time_of_message_clk", 
            val=np.take(self.dset_edit.time.gps.datetime, dset_has_idx),
            scale="gps",
            fmt="datetime",
            write_level="analysis",
        )
        dset.add_time(
            "has_time_of_message_clk", 
            val=np.take(self.dset_edit.tom.gps.datetime, dset_has_idx),
            scale="gps",
            fmt="datetime",
            write_level="analysis",
        )
        
        
    def add_code_bias_to_dataset(self, dset: "Dataset") -> None:
        """Add HAS code bias to dataset
        
        Args: 
            dset:   A Dataset containing model data.
        """
        # Relation between HAS signal definition and dataset naming, where RINEX convention are used
        #
        signal_def = {
            "E1-C": ("E","C1C"),
            "E5a-Q": ("E","C5Q"),
            "E6-C": ("E","C6C"), 
            "E5b-Q": ("E","C7Q"),
            "L1 C/A": ("G","C1C"),
            "L2 CL": ("G","C2L"),
            "L2 P": ("G","C2P"),
        }
    
        # Clean orbits by removing unavailable satellites, unhealthy satellites and checking validity length of
        # navigation records
        cleaners.apply_remover("gnss_clean_orbit_has", dset, orbit=self)

        # Get correct HAS message for given observations times and signal by determining the indices to HAS message
        # dataset
        #TODO-MURKS: This should be improved, that several GNSS can be handled by this function!!!
        if len(dset.unique("system")) > 1:
            log.fatal(f"Function add_code_bias_to_dataset() can not handle observations from two GNSS ({dset.unique('system')}).")

        for signal in set(self.dset_edit.signal):
            system, rinex_code = signal_def[signal]
            if system not in dset.unique("system"): # Skip not available GNSS
                continue
            if signal not in signal_def.keys():
                log.fatal(f"Galileo HAS code bias for signal {signal} handling is not implemented so far in Where.")

            dset_has_idx = self._get_has_message_idx(dset, system=system, signal=signal)
            dset.add_float(
                f"has_code_bias_{rinex_code.lower()}", 
                val=self.dset_edit.code_bias[dset_has_idx],
                unit="meters", 
                write_level="operational",
            )
        
        # Add additional information (TOD0: is that correct or can code bias information for each signal be different?)
        dset.add_float(
            "has_gnssiod_code_bias", 
            val=self.dset_edit.gnssiod[dset_has_idx], 
            unit="", 
            write_level="operational",
        )
        dset.add_float(
            "has_validity_code_bias", 
            val=self.dset_edit.validity[dset_has_idx],
            unit="",
            write_level="analysis",
        )
        dset.add_time(
            "has_reception_time_of_message_code_bias", 
            val=self.dset_edit.time.gps.datetime[dset_has_idx],
            scale="gps",
            fmt="datetime",
            write_level="analysis",
        )
        dset.add_time(
            "has_time_of_message_code_bias", 
            val=self.dset_edit.tom.gps.datetime[dset_has_idx],
            scale="gps",
            fmt="datetime",
            write_level="analysis",
        ) 


    def apply_code_bias_to_dataset(self, dset: "Dataset") -> None:
        """Apply HAS code bias to corresponding observation types in dataset
        
        Args: 
            dset:   A Dataset containing model data.
        """
        dset.meta.setdefault("has_corrected_obstypes", dict())
        
        # Clean orbits by removing unavailable satellites, unhealthy satellites and checking validity length of
        # navigation records
        cleaners.apply_remover("gnss_clean_orbit_has", dset, orbit=self)

        # Determine HAS code bias correction    
        for obs_type in sorted(dset.obs.fields):
            values = np.zeros(dset.num_obs)

            # Skip Doppler, carrier-phase and signal-to-noise observations types
            if obs_type.startswith("D") or obs_type.startswith("L") or obs_type.startswith("S"):
                continue

            for sys in dset.unique("system"):
                idx = dset.filter(system=sys)
  
                # Skip observation types, which are not defined for this GNSS
                if obs_type not in dset.meta["obstypes"][sys]:
                    continue

                sys_idx = self.dset_edit.filter(system=sys)
                available_signals = list(set(self.dset_edit.signal[sys_idx]))

                has_signal_type = self._get_has_bias_signal(sys, obs_type, available_signals)
                                
                if has_signal_type:
                    log.info(f"Use HAS code bias correction of signal '{sys}:{has_signal_type}' for " 
                             f"observation type '{sys}:{obs_type}'")
                    
                    # Get correct HAS message for given observations times by determining the indices to HAS message dataset
                    dset_has_idx = self._get_has_message_idx(dset, system=sys, signal=has_signal_type)
                    
                    values[idx] = self.dset_edit.code_bias[dset_has_idx][idx]

                    # Add HAS code bias correction information to meta variable
                    dset.meta["has_corrected_obstypes"].setdefault(sys, dict()).setdefault(has_signal_type, list())
                    dset.meta["has_corrected_obstypes"][sys][has_signal_type].append(obs_type)

                else:
                    log.warn(f"No Galileo HAS code bias correction are available for GNSS '{sys}' and observation "
                             f"type '{obs_type}'. Only for following signals are HAS code bias correction are "
                             f"available: {', '.join(available_signals)}")

            dset.add_float(
                f"has_code_bias.{obs_type}", 
                val=values, 
                unit="meters",
                write_level="analysis",
            )

            # Apply HAS code bias correction
            dset.obs[obs_type][:] = dset.obs[obs_type] + values
            

    def apply_phase_bias_to_dataset(self, dset: "Dataset") -> None:
        """Apply HAS phase bias to corresponding observation types in dataset
        
        Args: 
            dset:   A Dataset containing model data.
        """
        #TODO: Can functions apply_code_bias_to_dataset and apply_phase_bias_to_dataset function be merged together?

        dset.meta.setdefault("has_corrected_obstypes", dict())
        
        # Clean orbits by removing unavailable satellites, unhealthy satellites and checking validity length of
        # navigation records
        cleaners.apply_remover("gnss_clean_orbit_has", dset, orbit=self)

        # Determine HAS phase bias correction    
        for obs_type in dset.obs.fields:
            values = np.zeros(dset.num_obs)

            # Skip pseudo-range, Doppler and signal-to-noise observations types
            if obs_type.startswith("C") or obs_type.startswith("D") or obs_type.startswith("S"):
                continue

            for sys in dset.unique("system"):
                idx = dset.filter(system=sys)

                # Skip observation types, which are not defined for this GNSS
                if obs_type not in dset.meta["obstypes"][sys]:
                    continue

                sys_idx = self.dset_edit.filter(system=sys)
                available_signals = list(set(self.dset_edit.signal[sys_idx]))

                has_signal_type = self._get_has_bias_signal(sys, obs_type, available_signals)
                                
                if has_signal_type:
                    log.info(f"Use HAS phase bias correction of signal '{sys}:{has_signal_type}' for " 
                             f"observation type '{sys}:{obs_type}'")
                    
                    # Get correct HAS message for given observation times by determining the indices to HAS message dataset
                    dset_has_idx = self._get_has_message_idx(dset, system=sys, signal=has_signal_type)
                    
                    values[idx] = self.dset_edit.phase_bias[dset_has_idx][idx]

                    # Add HAS phase bias correction information to meta variable
                    dset.meta["has_corrected_obstypes"].setdefault(sys, dict()).setdefault(has_signal_type, list())
                    dset.meta["has_corrected_obstypes"][sys][has_signal_type].append(obs_type)

                else:
                    log.warn(f"No Galileo HAS phase bias correction are available for GNSS '{sys}' and observation "
                             f"type '{obs_type}'. Only for following signals are HAS phase bias correction are "
                             f"available: {', '.join(available_signals)}")

            dset.add_float(
                f"has_phase_bias.{obs_type}", 
                val=values, 
                unit="meters",
                write_level="analysis",
            )

            # Apply HAS phase bias correction
            dset.obs[obs_type][:] = dset.obs[obs_type] + values

    
    #
    # AUXILIARY FUNCTIONS
    #

    def _get_has_bias_signal(self, system: str, obs_type: str, available_signals: List[str]) -> Union[str, None]:
        """Get for given RINEX observation type and system the related HAS bias signal

        Args:
            system:            GNSS system identifier
            obs_type:          RINEX observation type
            available_signals: Available HAS signals

        Returns:
            HAS bias signal or if not given
        """
        # Relation between RINEX observation type and HAS signal definition
        #
        # Note: The HAS Galileo code bias corrections for signal E1-C, E5a-Q, E6-C, E5b-Q, L1 C/A, L2 CL and L2 P are 
        #       used for all defined RINEX observation types, which belongs to the same frequency.
        signal_def = {
            "E": {
                "E1-C": ["C1A", "C1B", "C1C", "C1X", "C1Z", 
                         "L1A", "L1B", "L1C", "L1X", "L1Z"],
                "E5a-Q": ["C5I", "C5Q", "C5X", "L5I", "L5Q", "L5X"],
                "E5b-Q": ["C7I", "C7Q", "C7X", "L7I", "L7Q", "L7X"],
                "E6-C": ["C6A", "C6B", "C6C", "C6X", "C6Z", 
                         "L6A", "L6B", "L6C", "L6X", "L6Z", ], 
                #"E??": ["C8I", "C8Q", "C8X", "L8I", "L8Q", "L8X"], # TODO: Not defined so ar.
            },
            "G": {
                "L1 C/A": ["C1C", "C1S", "C1L", "C1X", "C1P", "C1W", "C1Y", "C1M", 
                           "L1C", "L1S", "L1L", "L1X", "L1P", "L1W", "L1Y", "L1M", "L1N"],
                "L2 CL": ["C2S", "C2L", "C2X", "L2S", "L2L", "L2X"],
                "L2 P": ["C2C", "C2D", "C2S", "C2L", "C2X", "C2P", "C2W", "C2Y", "C2M", 
                         "L2C", "L2D", "L2S", "L2L", "L2X", "L2P", "L2W", "L2Y", "L2M", "L2N"],
            },
        }

        if system == "E":
            signals = ["E1-C", "E5a-Q", "E5b-Q", "E6-C"]
        elif system == "G":
            signals = ["L1 C/A", "L2 CL", "L2 P"]
        else:
            log.warn(f"For GNSS {system} are no Galileo HAS bias corrections provided.")
            return None

        for signal in signals:
            if signal in available_signals:
                if obs_type in signal_def[system][signal]:
                    return signal

        return None


    # TODO: Add functionality to .midgard/data/_position.py
    def _ric2trs(self, posvel: "PosVel") -> np.ndarray:
        """Transformation matrix from satellite coordinate system given with radial, in-track and cross-track
        directions to ITRS.
        
        Args: 
            posvel: Object of PosVel class
        
        Return:
            Transformation matrix
        """ 
        i_unit = posvel.trs.vel.unit_vector
        c_unit = nputil.unit_vector(np.cross(posvel.trs.pos, posvel.trs.vel))
        r_unit = np.cross(i_unit, c_unit)

        return np.stack((r_unit, i_unit, c_unit), axis=2)
    
    
    # TODO: Add functionality to .midgard/data/_position.py
    def _trs2ric(self, posvel: "PosVel") -> np.ndarray:
        """Transformation matrix from ITRS to satellite coordinate system given with radial, in-track and cross-track.
        
        Args: 
            posvel: Object of PosVel class
        
        Return:
            Transformation matrix
        """
        if posvel.ndim == 1:
            ric2trs = self._trs2ric(posvel).T
        else:
            ric2trs = self._trs2ric(posvel).transpose(0, 2, 1)
        return ric2trs


    @staticmethod
    def _add_dim(array: np.ndarray) -> np.ndarray:
        """Add dimension to 0 dimensional array
        
        If dimension of Numpy array is zero, then a dimension is added.
        
        Args:
            array: Numpy array with 0 or higher dimension
            
        Returns:
            Numpy array with at least 1 dimension
        """
        if type(array) == np.ndarray:  # Handling of numpy arrays
            output = np.expand_dims(array, axis=0) if array.ndim == 0 else array
        else: # Handling of scalar values
            output = np.expand_dims(array, axis=0)
            
        return output
    
    
    def _get_has_message_idx(
                self, 
                dset: "Dataset", 
                time: str = "time", 
                system: Union[str, None] = None,
                signal: Union[str, None] = None,
        ) -> List[int]:
        """Get Galileo HAS message indices for given observation epochs

        The indices relate the observation epoch to the correct set of HAS message. First the time difference
        between the observation epoch and a selected time is calculated to determine the correct HAS message. The 
        seleted time can be either the HAS message transmission time (receiver reception time) or the reference
        time of HAS message (TOM). Afterwards the HAS message is selected with the smallest time difference.
        Following option can be choosen for configuration file option 'has_message_nearest_to':

        | Option                       | Description                                                                  |
        |------------------------------|------------------------------------------------------------------------------|
        | tom                          | HAS message for given observation epoch is selected nearest to reference     |
        |                              | time of HAS message (TOM).                                                   |
        | tom:positive                 | Same as 'tom' option, but the difference between observation epoch and TOM   |
        |                              | has to be positive.                                                          |
        | transmission_time            | HAS message for given observation epoch is selected nearest to transmission  |
        |                              | time (receiver reception time).                                              |
        | transmission_time:positive   | Same as 'transmission_time' option, but the difference between observation   |
        |                              | epoch and transmission time has to be positive.                              |
        
        Satellite number and GNSS IOD is needed to select correct HAS message in relation to used broadcast navigation
        messages. In addition filtering related to given GNSS system identifier and signal is needed for code and 
        phase bias information.

        Args:
            dset:   A Dataset containing model data.
            time:   Define time fields to be used. It can be for example 'time' or 'sat_time'. 'time' is related to 
                    observation time and 'sat_time' to satellite transmission time.
            system: Satellite system used to get correct code or phase bias HAS messages in combination with 'signal' definition
            signal: Satellite signal used to get correct code or phase bias HAS messages in combination with 'system' definition

        Returns:
            HAS messages indices for given observation epochs.
        """
        has_message_nearest_to_options = [
            "tom",
            "tom:positive",
            "transmission_time",
            "transmission_time:positive",
        ]

        # Get configuration option
        has_message_nearest_to = config.tech.get("has_message_nearest_to", default="transmission_time:positive").str.rsplit(":", 1)
        if ":".join(has_message_nearest_to) not in has_message_nearest_to_options:
            log.fatal(
                f"Unknown value {':'.join(has_message_nearest_to)!r} for configuration option 'has_message_nearest_to'. "
                f"The following values can be selected: {', '.join(has_message_nearest_to_options)}"
            )

        time_key = has_message_nearest_to[0]
        positive = True if "positive" in has_message_nearest_to else False

        log.debug(f"HAS message is selected nearest to '{'+' if positive else '+/-'}{time_key}' time.")
        
        # Transmission time is 'time' field
        if time_key == "transmission_time":
            time_key = "time"

        # Check if HAS message are available
        not_available_sat = sorted(set(dset.satellite) - set(self.dset_edit.satellite))
        if not_available_sat:
            log.warn(
                f"The following satellites are not given in apriori HAS message file "
                f"{', '.join(self.dset_edit.meta['parser']['file_path'])}: {', '.join(not_available_sat)}"
            )

            cleaners.apply_remover("ignore_satellite", dset, satellites=not_available_sat)

        # Determine HAS message index for a given satellite, observation epoch and eventually system/signal
        indices = np.full(dset.num_obs, -1, dtype=int) # -1 is chosen to guarantee that not a wrong index is used 
                                                       # for getting correct HAS message

        has_epoch_mjd = self._add_dim(self.dset_edit[time_key].gps.mjd) 
        for idx_dset, (sat, obs_epoch) in enumerate(zip(dset.satellite, dset[time])):

            # Skip not relevant GNSS satellites, if GNSS system is defined
            if system:
                if not sat[0:1] == system:
                    continue
                        
            if signal:
                idx_has = self.dset_edit.filter(satellite=sat, system=system, signal=signal)
               
                if np.any(idx_has) == False:
                    # TODO-MURKS: Better solution for handling of L2 CL signal type. L2 CL bias corrections are not available for all satellites.
                    if system == "G" and signal == "L2 CL":
                        idx_has = self.dset_edit.filter(satellite=sat, system=system, signal="L2 P")

                        if np.any(idx_has) == False:
                            log.fatal(f"No valid HAS message could be found for satellite {sat}, GNSS signal {system}:{signal} and "
                                      f"observation epoch {obs_epoch.isot}. Use 'gnss_clean_orbit_has' remover.")

                        log.debug(f"HAS bias correction for signal 'G:{signal}' does not exists for satellite {sat}, "
                                  f"therefore use of corrections for signal 'G:L2 P' satellite {sat}.")

                    else:
                        log.fatal(f"No valid HAS message could be found for satellite {sat}, GNSS signal {system}:{signal} and "
                                  f"observation epoch {obs_epoch.isot}. Use 'gnss_clean_orbit_has' remover.")            
            else: 
                idx_has = self.dset_edit.filter(satellite=sat)
               
                if np.any(idx_has) == False:
                    log.fatal(f"No valid HAS message could be found for satellite {sat}, and observation epoch "
                             f"{obs_epoch.isot}. Use 'gnss_clean_orbit_has' remover.")

            nearest_idx = self._get_nearest_idx(has_epoch_mjd[idx_has], obs_epoch, time_key, positive) 
            indices[idx_dset] = idx_has.nonzero()[0][nearest_idx]

        return indices

    
    def _get_nearest_idx(self, has_epoch_mjd: np.ndarray, obs_epoch: "Time", time_key: str, positive: bool) -> np.ndarray: 
        """Get nearest HAS message data index for given observation epoch
        
        Args:
            has_epoch_mjd:  HAS message epoch filtered e.g. after satellite and GNSS IOD
            obs_epoch: 	    Observation epoch as Time object
            time_key:       Time key
            positive:       Difference between observation epoch and HAS message epoch has to be positive

        Returns:
            Nearest HAS messages indices for given observation epochs
            for given observation epoch
        """
        diff = obs_epoch.gps.mjd - has_epoch_mjd

        if positive:
            data = np.array([99999 if v < 0  else v for v in diff])
            if np.all(data == 99999): # No HAS message epochs larger than observation epoch
                log.fatal(f"No valid HAS message could be found before observation epoch {obs_epoch.gps.isot} (nearest HAS " 
                          f"message receiver reception time: {min(self.dset_edit[time_key].gps.isot)})")
            
            nearest_idx = data.argmin()
        else:
            nearest_idx = np.array([abs(diff)]).argmin()
            
        return nearest_idx
        

        


