"""Handling of GNSS broadcast orbits

Description:
------------
The module includes a class for handling apriori GNSS broadcast orbits.

Example:
    from where import apriori

    # Get broadcast ephemeris object
    brdc = apriori.get('orbit', rundate=rundate, station=station, , system=system, apriori_orbit='broadcast')

    # Write calculated Dataset to file
    brdc.dset.write()

"""

# Standard library imports
from datetime import timedelta
from typing import Dict, List, Union

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.collections import enums
from midgard.dev import plugins
from midgard.files import dependencies
from midgard.math.constant import constant
from midgard.math.unit import Unit
from midgard.parsers import rinex_nav

# Where imports
from where.data import dataset3 as dataset
from where import cleaners
from where.apriori import orbit
from where.lib import config
from where.lib import log
from where.lib import rotation

# Earth's gravitational constant
GM = {
    "C": constant.get("GM", source="cgcs2000"),
    "E": constant.get("GM", source="gtrf"),
    "G": constant.get("GM", source="wgs84"),
    "I": constant.get("GM", source="wgs84"),
    "J": constant.get("GM", source="jgs"),
}

# Earth's rotation rate
OMEGA = {
    "C": constant.get("omega", source="cgcs2000"),
    "E": constant.get("omega", source="gtrf"),
    "G": constant.get("omega", source="wgs84"),
    "I": constant.get("omega", source="wgs84"),
    "J": constant.get("omega", source="jgs"),
}


@plugins.register
class BroadcastOrbit(orbit.AprioriOrbit):
    """A class for representing apriori broadcast orbits

    RINEX navigation files can be read and edited. In addition GNSS satellite position, velocities, satellite clock
    correction and relativistic clock corrections can be calculated for given observation epochs based on broadcast
    ephemeris.

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
        relativistic_clock_correction():  Determine relativistic clock correction due to orbit eccentricity
        satellite_clock_correction():  Determine satellite clock correction
        unhealthy_satellites(): Get unhealthy satellites based on RINEX navigation file information
        _calculate():           Calculate broadcast ephemeris and satellite clock correction for given observation
                                epochs
        _edit():                Edit RINEX navigation file data and save it in a Dataset
        _read():                Read RINEX navigation file data and save it in a Dataset
    """

    name = "broadcast"

    def __init__(self, rundate, system, station=None, file_key=None, file_path=None, day_offset=1):
        """Set up a new BroadcastOrbit object, does not parse any data

        TODO: Remove dependency on rundate, use time to read correct files. (What to do with dataset?)

        Args:
            rundate (date):     Date of model run.
            station (str):      4 character station identifier.
            system (tuple):     List with GNSS system string codes (G, E, R, etc.).
            file_key (str):     Key to the broadcast orbit file defined in files.conf file.
            file_path (pathlib.PosixPath):  File path to broadcast orbit file.
            day_offset (int):   Day offset used to calculate the number of days to read.
        """
        super().__init__(rundate=rundate)
        self.system = system
        self.file_key = "gnss_rinex_nav_{system}" if file_key is None else file_key
        self.file_path = file_path
        self.day_offset = day_offset

        # TODO hjegei: Should it be not enough to 'station' in _dset_raw?
        if station:
            self._dset_raw.vars["station"] = station.lower()
            self._dset_raw.vars["STATION"] = self._dset_raw.vars["station"].upper()
            self._dset_edit.vars["station"] = station.lower()
            self._dset_edit.vars["STATION"] = self._dset_raw.vars["station"].upper()


    def _read(self, dset_raw):
        """Read RINEX navigation file data and save it in a Dataset

        Note that beside the given day also the navigation files from the day before and the day after is read.

        One RINEX navigation file is normally written for each GNSS. The navigation file extension depends on the GNSS
        (GPS: *.n, Galileo: *.l, GLONASS: *.g, ...). Therefore we have defined for each GNSS the navigation file
        name in the Where configuration file `files.conf`. In addition mixed navigation files exists, which includes
        navigation messages of different GNSS. We use following file keys in `files.conf`:
           
           | System   |  File key         |
           | :------- | :---------------- |
           | Galileo  | gnss_rinex_nav_E  |
           | GLONASS  | gnss_rinex_nav_R  |
           | GPS      | gnss_rinex_nav_G  |
           | Mixed    | gnss_rinex_nav_M  |

        Depending on the configuration options `systems` and `use_mixed_brdc_file` following navigation files are read:

         | Option                | File key          | What kind of navigation file is read? |
         | :---------------------| :-----------------| :-------------------------------------|
         | systems = G           | gnss_rinex_nav_G  | Only the GPS navigation file          |
         | systems = G E         | gnss_rinex_nav_G  | GPS and Galileo navigation files      |
         |                       | gnss_rinex_nav_E  |                                       |
         | use_mixed_brdc_file   | gnss_rinex_nav_M  | Mixed GNSS navigation file            |

        Args:
            dset_raw (Dataset):     Dataset representing raw data from RINEX navigation file
        """
        date_to_read = dset_raw.analysis["rundate"] - timedelta(days=self.day_offset)
        use_mixed_brdc_file = config.tech.get("use_mixed_brdc_file", default=False).bool
        systems = {"M"} if use_mixed_brdc_file == True else self.system
        file_paths = list()

        # Loop over days to read
        while date_to_read <= dset_raw.analysis["rundate"] + timedelta(days=self.day_offset):

            date = date_to_read.strftime("%Y-%m-%d")
            meta = dict()

            for sys in systems:
                if self.file_path is None:
                    file_path = config.files.path(
                        self.file_key.format(system=sys), file_vars=config.date_vars(date_to_read)
                    )
                else:
                    file_path = self.file_path

                log.debug(f"Parse broadcast orbit file {file_path}")

                # Generate temporary Dataset with orbit file data
                dset_temp = dataset.Dataset(
                    rundate=date_to_read, pipeline=dset_raw.vars["pipeline"], stage="temporary"
                )
                parser = rinex_nav.get_rinex2_or_rinex3(file_path)
                dset_temp.update_from(parser.as_dataset())
                file_paths.append(str(parser.file_path))
                dependencies.add(str(parser.file_path), label=self.file_key.format(system=sys))  # Used for output writing

                # Extend Dataset dset_raw with temporary Dataset
                if dset_raw.num_obs:

                    # Merge meta data information
                    # TODO: Handle meta data information correctly. Meta data information based on different GNSS navigation
                    #       message files has to be merged correctly together. What to do with 'sat_sys' and 'leap_seconds'?
                    if date in dict(dset_raw.meta).keys():
                        for key in ["iono_para", "time_sys_corr"]:
                            dset_temp.meta[key].update(dset_raw.meta[date][key])
                    dset_raw.extend(dset_temp, meta_key=date)

                else:
                    # Add date key to meta TODO: Better solution?
                    meta.setdefault(date, dict()).update(dset_temp.meta)
                    for k in dict(dset_temp["meta"]).keys():
                        del dset_temp.meta[k]
                    dset_temp.meta.update(meta)

                    dset_raw.update_from(dset_temp)

            date_to_read += timedelta(days=1)

        dset_raw.meta.add("file_path", file_paths, section="parser")

        if "E" in dset_raw.unique("system"):
            self._galileo_signal_health_status()


    def _edit(self, dset_edit: "Dataset") -> "Dataset":
        """Edit RINEX navigation file data and save it in a Dataset

        First the navigation messages are sorted after the satellite and time of transmission. Afterwards duplicated
        navigation messages for a satellite are removed, whereby the first occurrence of the navigation message is
        kept. The satellite navigation epoch (time of clock (toc)) and the IODE navigation number is used for
        indication of a duplicated epoch. The meaning of the IODE navigation number depends on the GNSS:

            - GPS: Ephemeris issue of data (IODE)
            - Galileo: Issue of Data of the NAV batch (IODnav)
            - QZSS: Ephemeris issue of data (IODE)
            - BeiDou: Age of Data Ephemeris (AODE)
            - IRNSS: Issue of Data, Ephemeris and Clock (IODEC)

        It can be defined with the configuration option 'navigation_message_type', what kind of navigation message type
        should be used in the analysis. Other navigation message types are removed from the broadcast ephemeris. For
        example for Galileo INAV or FNAV navigation messages can be chosen.

        Args:
            dset_edit (Dataset):     Dataset representing edited data from RINEX navigation file
        """

        # TODO: This does not work. -> Can not find remover ignore_duplicated_navigation_messages().
        # cleaners.removers.ignore_duplicated_navigation_messages(dset_edit) #MURKS

        # Generate Pandas frame for all navigation message entries
        # TODO: self.dset_raw.as_dataframe() could also be used, but performance is quite bad.
        nav = pd.DataFrame()
        for field in self.dset_raw.fields:
            if "midgard.data._time" in str(type(self.dset_raw[field])): #Handling of "Time" objects
                nav[field] = self.dset_raw[field].gps.iso
            nav[field] = self.dset_raw[field] 

        # Filter frame
        nav_filtered = nav.sort_values(by=["satellite", "time", "transmission_time"])

        # Remove duplicated navigation message entries (IODEs)
        # TODO: Maybe it should be a configuration option, how to filter duplicated epochs. Keep first and last.
        idx = nav_filtered.duplicated(
            subset=["satellite", "time", "iode", "nav_type", "transmission_time"], keep="first"
        )
        nav_duplicates = nav_filtered[["satellite", "time", "iode", "nav_type"]][idx]
        with pd.option_context("display.max_rows", None, "display.max_columns", 5):
            log.info(f"Remove {len(nav_duplicates)} duplicated navigation message entries.")
            log.debug(f"List of duplicated navigation messages: \n{nav_duplicates}")

        nav_filtered = nav_filtered.drop_duplicates(
            subset=["satellite", "time", "iode", "nav_type", "transmission_time"], keep="first"
        )

        # Remove navigation message types, which are not needed
        nav_type = config.tech.get("navigation_message_type", default="").dict
        if nav_type:
            for sys, type_ in nav_type.items():
                sys = sys[0] if len(sys) > 0 else sys

                # Overwrite message type definition if only 'INAV' or 'FNAV' is specified
                if type_ == "INAV":
                    type_ = ["INAV_E1", "INAV_E5b", "INAV_E1E5b"]
                elif type_ == "FNAV":
                    type_ = ["FNAV_E5a"]

                remove_nav_type = nav_filtered.query("system == @sys and nav_type != @type_")
                if not remove_nav_type.empty:
                    log.info(
                        f"Remove {', '.join(set(remove_nav_type.nav_type))!r} navigation messages of GNSS {sys!r}"
                    )
                    nav_filtered = pd.concat([nav_filtered, remove_nav_type]).drop_duplicates(keep=False)

        if nav_filtered.empty:
            log.fatal(f"No navigation messages available for GNSS: {','.join(self.dset_raw.unique('system'))} ")

        # TODO hjegei: Possible future ...
        # dset_edit.copy_from(self.dset_raw)
        # dset_edit.reorder(nav_filtered.index.values)

        # Convert edited fields to Dataset
        nav_np = nav_filtered.values
        fields = nav_filtered.columns

        dset_edit.vars["orbit"] = self.name
        dset_edit.num_obs = nav_filtered.shape[0]
        dset_edit.meta.update(self.dset_raw.meta)

        for idx, field in enumerate(fields):
            if field.startswith("time_") or field.startswith("transmission_time_") or field.startswith("toe_"):
                continue # Skip unnecessary fields
            elif field in ["time", "transmission_time", "toe"]:
                dset_edit.add_time(field, val=nav_np[:, idx], scale="gps", fmt="datetime")
            elif field in ["nav_type", "satellite", "system"]:
                dset_edit.add_text(field, val=nav_np[:, idx])
            else:
                unit = self.dset_raw.unit(field)[0] if self.dset_raw.unit(field) else None
                dset_edit.add_float(field, val=nav_np[:, idx].astype(float), unit=unit)


    def _calculate(self, dset_out: "Dataset", dset_in: "Dataset", time: str = "time") -> None:
        """Calculate broadcast ephemeris and satellite clock correction for given observation epochs

        As a first step observations are removed from unavailable satellites, unhealthy satellites and for exceeding the 
        validity length of navigation records. The input Dataset contains observation epochs for which the broadcast
        ephemeris and satellite clock correction should be determined. 

        Args:
            dset_out: Output Dataset representing calculated broadcast ephemeris with following fields:

        | Field                   | Type          | Unit    | Description                                             |
        | :-----------------------| :------------ | :------ | :------------------------------------------------------ |
        | gnss_satellite_clock    | numpy.ndarray | m       | Satellite clock correction                              |
        | gnss_relativistic_clock | numpy.ndarray | m       | Relativistic clock correction due to orbit eccentricity |
        | sat_posvel              | PosVelTable   | m       | Satellite position and velocity                         |
        | satellite               | numpy.ndarray |         | Satellite numbers                                       |
        | system                  | numpy.ndarray |         | GNSS identifiers                                        |
        | time                    | TimeTable     |         | Observation epochs                                      |
        | used_iode               | numpy.ndarray |         | IODE of selected broadcast ephemeris block              |
        | used_transmission_time  | TimeTable     |         | Transmission time of selected broadcast ephemeris block |
        | used_toe                | TimeTable     |         | Time of ephemeris (TOE) of selected broadcast ephemeris |
        |                         |               |         | block                                                   |
        
            dset_in:  Input Dataset containing model data for which broadcast ephemeris should be determined.
            time:     Define time fields to be used. It can be for example 'time' or 'sat_time'. 'time' is related to 
                      observation time and 'sat_time' to satellite transmission time.
        """
        not_implemented_sys = set(dset_in.system) - set("EG")
        if not_implemented_sys:
            log.warn(
                f"At the moment Where can provide broadcast ephemeris for GNSS 'E' and 'G', "
                f"but not for {', '.join(not_implemented_sys)}."
            )

            cleaners.apply_remover("gnss_ignore_system", dset_in, systems=not_implemented_sys)

        log.info(
            f"Calculating satellite position/velocity (broadcast) based on RINEX navigation file "
            f"{', '.join(self.dset_edit.meta['parser']['file_path'])}"
        )

        apply_has_correction = config.tech.get("apply_has_correction", default=False).bool
        
        # Get navigation block index for given observation epochs and clean orbits by removing unavailable satellites, 
        # unhealthy satellites and checking validity length of navigation records
        if "navigation_idx" not in dset_in.fields:
            cleaners.apply_remover("gnss_clean_orbit", dset_in, orbit_flag="broadcast", orbit=self)
        dset_brdc_idx = dset_in.navigation_idx.astype(int)

        # Add used navigation message type (useful for postprocessing of Where results)
        navigation_message_type = config.tech.navigation_message_type.dict
        systems = dset_in.unique("system")
        for sys in list(navigation_message_type.keys()):
            if sys not in dset_in.unique("system"):
                del navigation_message_type[sys]  # Remove unused navigation message types
        dset_in.meta["navigation_message_type"] = navigation_message_type

        # Loop over all observations
        # TODO: Generation of vectorized solution, if possible?

        # BUG: Use of GPSSEC does not work for GPS WEEK crossovers. MJD * Unit.day2second() would a better solution.
        #      The problem is that use of GPSSEC compared to MJD * Unit.day2second() is not consistent!!!!
        sat_pos = np.zeros((dset_in.num_obs, 3))
        sat_vel = np.zeros((dset_in.num_obs, 3))
        for obs_idx, (time_gpsweek, time_gpssec, brdc_idx, sys) in enumerate(
            # zip(dset_in[time].gps.gpsweek, dset_in[time].gps.gpssec, dset.navigation_idx.as_type(int), dset_in.system)
            zip(dset_in[time].gps.gps_ws.week, dset_in[time].gps.gps_ws.seconds, dset_brdc_idx, dset_in.system)
        ):

            # TODO: get_row() function needed for brdc -> brdc.get_row(kk)
            sat_pos[obs_idx], sat_vel[obs_idx] = self._get_satellite_position_velocity(
                time_gpsweek, time_gpssec, brdc_idx, sys
            )

            # +DEBUG
            # print("DEBUG: {} obs_idx: {:>5d} brdc_idx: {:>5d} toc: {:>5.0f} {:>6.0f} toe: {:>6.0f} trans_time: {:>6.0f}"
            #    " tk: {:>16.10f} iode: {:>3d} sqrt_a: {:>17.10f} sat_pos: {:>21.10f} {:>21.10f} {:>21.10f} "
            #    "sat_vel: {:>17.10f} {:>17.10f} {:>17.10f} sat_clk_bias: {:>17.10f}, sat_clk_drft: {:>17.10f} "
            #    ''.format(self.dset_edit.satellite[brdc_idx], obs_idx, brdc_idx,
            #    dset_in[time].gps.jd_frac[obs_idx] * 86400,
            #    dset_in[time].gps.gps_ws.seconds[obs_idx],
            #    self.dset_edit.toe.gps.gps_ws.seconds[brdc_idx],
            #    self.dset_edit.transmission_time.gps.gps_ws.seconds[brdc_idx],
            #    dset_in[time].gps.jd_frac[obs_idx]-self.dset_edit.toe.gps.gps_ws.seconds[brdc_idx],
            #    int(self.dset_edit.iode[brdc_idx]),
            #    self.dset_edit.sqrt_a[brdc_idx],
            #    sat_pos[obs_idx][0], sat_pos[obs_idx][1], sat_pos[obs_idx][2],
            #    sat_vel[obs_idx][0], sat_vel[obs_idx][1], sat_vel[obs_idx][2],
            #    self.dset_edit.sat_clock_bias[brdc_idx],
            #    self.dset_edit.sat_clock_drift[brdc_idx],)
            # )
            # -DEBUG

        # Copy fields from model data Dataset
        dset_out.num_obs = dset_in.num_obs
        dset_out.add_text("satellite", val=dset_in.satellite)
        dset_out.add_text("system", val=dset_in.system)
        dset_out.add_time("time", val=dset_in[time])
        dset_out.vars["orbit"] = self.name
        
        if apply_has_correction:
            dset_out.add_float("has_gnssiod_orb", val=dset_in.has_gnssiod_orb)

        # Add time field
        dset_out.add_time("used_transmission_time", val=self.dset_edit.transmission_time[dset_brdc_idx])
        dset_out.add_time("used_toe", val=self.dset_edit.toe[dset_brdc_idx])

        # Add float fields
        for field in ["bgd_e1_e5a", "bgd_e1_e5b", "tgd", "tgd_b1_b3", "tgd_b2_b3"]:
            if field in self.dset_edit.fields:
                dset_out.add_float(field, val=self.dset_edit[field][dset_brdc_idx])

        dset_out.add_float(
            "gnss_relativistic_clock", val=self.relativistic_clock_correction(sat_pos, sat_vel), unit="meter"
        )
        dset_out.add_float(
            "gnss_satellite_clock", val=self.satellite_clock_correction(dset_in, time=time), unit="meter"
        )
        dset_out.add_float("used_iode", val=self.dset_edit.iode[dset_brdc_idx])

        # Add satellite position and velocity to Dataset
        dset_out.add_posvel("sat_posvel", time=dset_out.time, system="trs", val=np.hstack((sat_pos, sat_vel)))

        # +DEBUG
        # for num_obs  in range(0, dset_out.num_obs):
        #    print('DEBUG: ', dset_out.satellite[num_obs],
        #                     dset_out.time.gps.datetime[num_obs],
        #                     dset_out.time.gps.mjd_frac[num_obs]*24*3600,
        #                     dset_out.gnss_satellite_clock[num_obs],
        #                     dset_out.gnss_relativistic_clock[num_obs])
        # -DEBUG
    


    #
    # CLASS METHODS
    #  
    def add_satellite_clock_parameter(self, dset: "Dataset", time: str = "time") -> np.ndarray:
        """Add satellite clock parameters to dataset

        Args:
            dset: A Dataset containing model data.
            time: Define time fields to be used. It can be for example 'time' or 'sat_time'. 'time' is related to 
                  observation time and 'sat_time' to satellite transmission time.

        """
        if "navigation_idx" not in dset.fields:
            log.fatal("Dataset field 'navigation_idx' does not exists. Remover 'gnss_clean_orbit' should be used.")
        dset_brdc_idx = dset.navigation_idx.astype(int)

        # Add satellite clock parameters to dataset
        dset.add_float("sat_clock_bias", val=self.dset_edit.sat_clock_bias[dset_brdc_idx], write_level="analysis")
        dset.add_float("sat_clock_drift", val=self.dset_edit.sat_clock_drift[dset_brdc_idx], write_level="analysis")
        dset.add_float(
            "sat_clock_drift_rate", 
            val=self.dset_edit.sat_clock_drift_rate[dset_brdc_idx], 
            write_level="analysis",
        )

         
    def get_ephemeris_field(self, field, dset: "Dataset", time: str = "time") -> np.ndarray:
        """Get given ephemeris field values for belonging navigation records

        Args:
            field:  Field name of navigation record.
            dset:   A Dataset containing model data.
            time:   Define time fields to be used. It can be for example 'time' or 'sat_time'. 'time' is related to 
                    observation time and 'sat_time' to satellite transmission time.

        Returns:
            Values for given field
        """
        return self.dset_edit[field][dset.navigation_idx.astype(int)]


    def get_nearest_idx(
                self, 
                idx: np.ndarray, 
                obs_epoch: float, 
                time_key: str, 
                positive: bool,
                satellite: str,
                iode: Union[float, None] = None,
                log_level: str = "fatal",
        ) -> float:
            """Get nearest navigation message data index for given observation epoch
            
            Args:
                idx:        Index used to filter navigation message e.g. after satellite, GNSS IOD
                obs_epoch:  Observation epoch in modified Julian day
                time_key:   Time key
                positive:   Difference between observation epoch and broadcast navigation message has to be positive
                satellite:  Satellite name
                iode:       IODE, which should be used to select broadcast navigation message (e.g. relevant in case
                            to match GNSS IOD of Galileo HAS message against broadcast navigation message IODE)
                log_level:  Define log level in case no broadcast navigation message data are available 

            Returns:
                Nearest broadcast navigation messages indices for given observation epochs or None in case no 
                corresponding broadcast navigation message could be found (this is only valid if 'log_level=fatal')
            """
            # Select nearest navigation message related to a given IODE and observation epoch
            if iode:
                idx_iode = self.dset_edit.iode[idx] == iode
    
                diff = np.full(len(self.dset_edit.iode[idx]), 99999.0)
                diff[idx_iode] = obs_epoch - self._add_dim(self.dset_edit[time_key].gps.mjd)[idx][idx_iode]

            # Select nearest naviation message related to given observation epoch
            else: 
                diff = obs_epoch - self._add_dim(self.dset_edit[time_key].gps.mjd)[idx]
                
                
            if positive:
                data = np.array([99999 if v < 0 else v for v in diff])
                if np.all(data == 99999): # No broadcast navigation message epochs larger than observation epoch
                    getattr(log, log_level)(f"No valid broadcast navigation message could be found for satellite " 
                              f"{satellite} before observation epoch {obs_epoch} (nearest broadcast navigation message "
                              f"receiver reception time: {min(self.dset_edit[time_key].gps.mjd[idx])})") 
                    return None
                
                nearest_idx = data.argmin()
            else:
                nearest_idx = np.array([abs(diff)]).argmin()
                
            return nearest_idx 


    def satellite_clock_correction(self, dset: "Dataset", time: str = "time") -> np.ndarray:
        """Determine satellite clock correction

        The satellite clock correction is based on Section 20.3.3.3.3.1 in :cite:`is-gps-200h`.

        Args:
            dset: A Dataset containing model data.
            time: Define time fields to be used. It can be for example 'time' or 'sat_time'. 'time' is related to 
                  observation time and 'sat_time' to satellite transmission time.

        Returns:
            GNSS satellite clock corrections for each observation in [m] (Note: without relativistic orbit eccentricity
            correction)
        """
        # Get correct navigation block for given observations times by determining the indices to broadcast ephemeris
        # Dataset
        if "navigation_idx" not in dset.fields:
            dset_brdc_idx = self._get_brdc_block_idx(dset)
        else:
            dset_brdc_idx = dset.navigation_idx.astype(int)

        # Elapsed time referred to clock data reference epoch toc in [s]
        # BUG: Use of GPSSEC does not work for GPS WEEK crossovers. MJD * Unit.day2second() would a better solution. The
        #     problem is that use of GPSSEC compared to MJD * Unit.day2second() is not consistent!!!!
        gpsweek_diff = (
            dset[time].gps.gps_ws.week - self._add_dim(self.dset_edit.time.gps.gps_ws.week)[dset_brdc_idx]
        ) * Unit.week2second
        tk = dset[time].gps.gps_ws.seconds - self._add_dim(self.dset_edit.time.gps.gps_ws.seconds)[dset_brdc_idx] + gpsweek_diff

        return (
            self.dset_edit.sat_clock_bias[dset_brdc_idx]
            + self.dset_edit.sat_clock_drift[dset_brdc_idx] * tk
            + self.dset_edit.sat_clock_drift_rate[dset_brdc_idx] * tk ** 2
        ) * constant.c
        
    def satellite_clock_rate_correction(self, dset: "Dataset", time: str = "time") -> np.ndarray:
        """Determine satellite clock rate correction

        The satellite clock rate correction is based on Section 6.2.6.2 in :cite:`zhang2007`.

        Args:
            dset: A Dataset containing model data.
            time: Define time fields to be used. It can be for example 'time' or 'sat_time'. 'time' is related to 
                  observation time and 'sat_time' to satellite transmission time.

        Returns:
            GNSS satellite clock rate corrections for each observation in [m]
        """
        # Get correct navigation block for given observations times by determining the indices to broadcast ephemeris
        # Dataset
        if "navigation_idx" not in dset.fields:
            dset_brdc_idx = self._get_brdc_block_idx(dset)
        else:
            dset_brdc_idx = dset.navigation_idx.astype(int)

        # Elapsed time referred to clock data reference epoch toc in [s]
        # BUG: Use of GPSSEC does not work for GPS WEEK crossovers. MJD * Unit.day2second() would a better solution. The
        #     problem is that use of GPSSEC compared to MJD * Unit.day2second() is not consistent!!!!
        gpsweek_diff = (
            dset[time].gps.gps_ws.week - self.dset_edit.time.gps.gps_ws.week[dset_brdc_idx]
        ) * Unit.week2second
        tk = dset[time].gps.gps_ws.seconds - self.dset_edit.time.gps.gps_ws.seconds[dset_brdc_idx] + gpsweek_diff

        return (
            self.dset_edit.sat_clock_drift[dset_brdc_idx]
            + 2 * self.dset_edit.sat_clock_drift_rate[dset_brdc_idx] * tk
        ) * constant.c


    def signal_health_status(
                self, 
                dset: "Dataset", 
                frequencies: Union[None, List] = None,
    ) -> np.ndarray:
        """Determine signal health status

        How the signal health status has to be handled depends on GNSS, used observation type and navigation message
        type. The determination of the signal health status is defined for:

        | System   | Type  | Section       | Reference                          |
        | :------- | :---- | :------------ | :--------------------------------- |
        | GPS      | LNAV  | 20.3.3.3.1.4  | :cite:`is-gps-200h`                |
        |          | CNAV  | 30.3.3.1.1.2  | :cite:`is-gps-200h`                |
        | Galileo  | all   | 2.3           | :cite:`galileo-os-sdd`             |

        
        The signal health status is indicated by different numbers:

        | Number | Status          | Description                                                                     |
        | :----- | :-------------- | :------------------------------------------------------------------------------ |
        | 0      | Healthy         |                                                                                 |
        | 1      | Marginal (SISA) | Defined only for Galileo and describes marginal status based on signal-in-space |
        |        |                 | accuracy (SISA)                                                                 |
        | 2      | Marginal (DVS)  | Defined only for Galileo and describes marginal status based on data validity   |
        |        |                 | status.                                                                         |
        | 3      | Unhealthy       |                                                                                 |

        TODO: Add handling of BeiDou, QZSS and IRNSS signal health status.

        Args:
            dset:           A Dataset containing model data.
            frequencies:    TODO: Why was this arguments introduced???

        Returns:
            Array with GNSS signal health status for each observation
        """
        signal_health_status = np.zeros(dset.num_obs, dtype=bool)
        nav_type = config.tech.get("navigation_message_type", default="").dict

        # Get correct navigation block for given observations times by determining the indices to broadcast ephemeris
        # Dataset
        if "navigation_idx" not in dset.fields:
            dset_brdc_idx = self._get_brdc_block_idx(dset)
        else:
            dset_brdc_idx = dset.navigation_idx.astype(int)

        for sys in dset.unique("system"):
            idx = dset.filter(system=sys)
            data_tmp = np.zeros(len(dset.system[idx]), dtype=bool)

            # TODO: If several editors are used can it happen, that GNSS is not given in dset.meta["obstypes"], but
            #      observations still existing.
            if "obstypes" in dset.meta.keys():
                if sys not in dset.meta["obstypes"].keys():
                    continue

            if sys in ["C", "G", "I", "J", "R"]:
                idx_obs = self.dset_edit.sv_health[dset_brdc_idx][idx] > 0
                data_tmp[idx_obs] = 3

            elif sys == "E":

                # Get frequencies by keeping a selection of E1 (1), E5a (5) and E5b (7) frequencies
                # TODO: It has to be checked, if following selection of frequencies is correct?
                if frequencies is None:
                    if nav_type["E"] == "FNAV":
                        frequencies = ["E5a"]
                    elif nav_type["E"] == "INAV":
                        keep_freq = {"1", "7"}

                        if dset.vars["pipeline"] == "sisre":
                            frequencies = config.tech.frequencies.dict["E"].split("_")
                        elif dset.vars["pipeline"] == "gnss":
                            freq_numbers = (
                                set([t[1] for t in dset.meta["obstypes"]["E"]]) & keep_freq
                            )  # 1: E1, 5: E5a, 7: E5b
                            frequencies = [enums.gnss_num2freq_E["f" + num] for num in freq_numbers]
                        else:
                            frequencies = ["E1", "E5b"]

                for freq in frequencies:
                    data_tmp = self.dset_edit["sis_status_" + freq.lower()][dset_brdc_idx][idx]

            else:
                log.fatal(f"GNSS '{sys}' is not defined.")
            signal_health_status[idx] = data_tmp

        # +DEBUG
        ## Add signal health status to dataset
        # if "signal_health_status" in dset.fields:
        #    dset.signal_health_status[:] = signal_health_status
        # else:
        #    dset.add_float("signal_health_status", val=signal_health_status, write_level="analysis"))
        # -DEBUG

        return signal_health_status


    def total_group_delay(self, dset: "Dataset") -> np.ndarray:
        """Determine total group delay

        What kind of total group delay has to be applied depends on GNSS, used observation type and navigation message
        type. The determination of the total group delay is defined for GPS in Section 20.3.3.3.3.2 in 
        :cite:`is-gps-200h` and for Galileo in Section 5.1.5 in :cite:`galileo-os-sis-icd`. 

        TODO: Add handling of BeiDou, QZSS and IRNSS TGDs.

        Args:
            dset: A Dataset containing model data.

        Returns:
            GNSS total group delay for each observation in [m] 
        """
        # Get correct navigation block for given observations times by determining the indices to broadcast ephemeris
        # Dataset
        if "navigation_idx" not in dset.fields:
            dset_brdc_idx = self._get_brdc_block_idx(dset)
        else:
            dset_brdc_idx = dset.navigation_idx.astype(int)

        for sys, obstype in dset.meta["obstypes"].items():

            idx = dset.filter(system=sys)

            if len(obstype) != 1:
                log.warn("Total group delay can only be applied for single-frequency solutions.")
                return np.zeros(dset.num_obs)

            freq = obstype[0][1]

            if sys == "G":
                f_L1 = enums.gnss_freq_G.L1
                f_L2 = enums.gnss_freq_G.L2

                if freq == "1":
                    tgd = self.dset_edit.tgd[dset_brdc_idx][idx] * constant.c
                elif freq == "2":
                    tgd = f_L1 ** 2 / f_L2 ** 2 * self.dset_edit.tgd[dset_brdc_idx][idx] * constant.c

            elif sys == "E":
                nav_type = config.tech.navigation_message_type.dict["E"]

                f_E1 = enums.gnss_freq_E.E1
                f_E5a = enums.gnss_freq_E.E5a
                f_E5b = enums.gnss_freq_E.E5b

                if nav_type.startswith("FNAV"):
                    if freq == "1":
                        tgd = self.dset_edit.bgd_e1_e5a[dset_brdc_idx][idx] * constant.c
                    elif freq == "5":
                        tgd = f_E1 ** 2 / f_E5a ** 2 * self.dset_edit.bgd_e1_e5a[dset_brdc_idx][idx] * constant.c

                elif nav_type.startswith("INAV"):
                    if freq == "1":
                        tgd = self.dset_edit.bgd_e1_e5b[dset_brdc_idx][idx] * constant.c
                    elif freq == "7":
                        tgd = f_E1 ** 2 / f_E5b ** 2 * self.dset_edit.bgd_e1_e5b[dset_brdc_idx][idx] * constant.c

        return tgd
   
     
    #
    # AUXILIARY FUNCTIONS
    #

    def _add_dim(self, array: np.ndarray) -> np.ndarray:
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


    def _galileo_signal_health_status(self):
        """Determine Galileo signal health status and add sis_status_<signal> fields to Dataset
       
        The determination of the signal health status is defined for Galileo in Section 2.3 of :cite:`galileo-os-sdd`.

        The signal health status is indicated by different numbers:

        | Number | Status          | Description                                               |
        |--------|-----------------|-----------------------------------------------------------|
        | 0      | Healthy         |                                                           |
        | 1      | Marginal (SISA) | Marginal status based on signal-in-space accuracy (SISA). |
        | 2      | Marginal (DVS)  | Marginal status based on data validity status.            |
        | 3      | Unhealthy       |                                                           |

        Args:
            dset: A Dataset containing model data.
        """

        data = dict()
        idx = self.dset_raw.filter(system="E")

        for signal in ["e1", "e5a", "e5b"]:

            field = "sis_status_" + signal
            data = np.zeros(self.dset_raw.num_obs)
            data_tmp = np.zeros(len(self.dset_raw.system[idx]))

            # Check Signal-in-space Accuracy (SISA)
            idx_obs = np.logical_or(self.dset_raw.sv_accuracy[idx] < 0, self.dset_raw.sv_accuracy[idx] >= 126)
            data_tmp[idx_obs] = 1

            # Check Data Validity Status (DVS)
            idx_obs = self.dset_raw["dvs_" + signal][idx] == True
            data_tmp[idx_obs] = 2

            # Check Signal Health Status (SHS)
            idx_obs = self.dset_raw["shs_" + signal][idx] == True
            data_tmp[idx_obs] = 3

            data[idx] = data_tmp
            self.dset_raw.add_float(field, val=data)


    def _get_brdc_block_idx(self, dset: "Dataset", time: str = "time") -> List[int]:
        """Get GNSS broadcast ephemeris block indices for given observation epochs

        The indices relate the observation epoch to the correct set of broadcast ephemeris. First the time difference
        between the observation epoch and a selected time is calculated to determine the correct broadcast ephemeris
        block. The seleted time can be either the navigation epoch (time of clock (TOC)), the time of ephemeris (TOE)
        or the transmission time. Afterwards the broastcast block is selected with the smallest time difference.
        Following option can be choosen for configuration file option 'brdc_block_nearest_to':

        | Option                       | Description                                                                  |
        | :--------------------------- | :--------------------------------------------------------------------------- |
        | toc                          | Broadcast block for given observation epoch is selected nearest to           |
        |                              | navigation epoch (time of clock (TOC)).                                      | 
        | toc:positive                 | Same as 'toc' option, but the difference between observation epoch and TOC   |
        |                              | has to be positive.                                                          |
        | toe                          | Broadcast block for given observation epoch is selected nearest to time of   |
        |                              | ephemeris (TOE).                                                             |
        | toe:positive                 | Same as 'toe' option, but the difference between observation epoch and TOE   |
        |                              | has to be positive.                                                          |
        | transmission_time            | Broadcast block for given observation epoch is selected nearest to           |
        |                              | transmission time.                                                           |
        | transmission_time:positive   | Same as 'transmission_time' option, but the difference between observation   |
        |                              | epoch and transmission time has to be positive.                              |

        Args:
            dset: A Dataset containing model data.
            time: Define time fields to be used. It can be for example 'time' or 'sat_time'. 'time' is related to 
                  observation time and 'sat_time' to satellite transmission time.

        Returns:
            Broadcast ephemeris block indices for given observation epochs.
        """
        brdc_block_nearest_to_options = [
            "toc",
            "toc:positive",
            "toe",
            "toe:positive",
            "transmission_time",
            "transmission_time:positive",
        ]
        brdc_idx = list()
        apply_has_correction = config.tech.get("apply_has_correction", default=False).bool

        # Get configuration option
        brdc_block_nearest_to = config.tech.get("brdc_block_nearest_to", default="toe:positive").str.rsplit(":", 1)
        if ":".join(brdc_block_nearest_to) not in brdc_block_nearest_to_options:
            log.fatal(
                f"Unknown value {':'.join(brdc_block_nearest_to)!r} for configuration option 'brdc_block_nearest_to'. "
                f"The following values can be selected: {', '.join(brdc_block_nearest_to_options)}"
            )

        time_key = brdc_block_nearest_to[0]
        positive = True if "positive" in brdc_block_nearest_to else False

        log.debug(f"Broadcast block is selected nearest to '{'+' if positive else '+/-'}{time_key}' time.")
        
        # Time of clock (TOC) is 'time' field
        if time_key == "toc":
            time_key = "time"

        # Check if broadcast orbits are available
        not_available_sat = sorted(set(dset.satellite) - set(self.dset_edit.satellite))
        if not_available_sat:
            log.warn(
                f"The following satellites are not given in apriori broadcast orbit file "
                f"{', '.join(self.dset_edit.meta['parser']['file_path'])}: {', '.join(not_available_sat)}"
            )

            cleaners.apply_remover("ignore_satellite", dset, satellites=not_available_sat)
                     
        if apply_has_correction:
            
            # Determine broadcast ephemeris block index for a given satellite, HAS IOD and observation epoch
            for sat, iode, obs_epoch in zip(dset.satellite, dset.has_gnssiod_orb, dset[time]):  # TODO: Is _add_time needed for dset[time]?
    
                idx = self.dset_edit.filter(satellite=sat, iode=iode)
                
                if np.any(idx) == False:
                    log.fatal(f"No valid broadcast navigation message could be found for satellite {sat}, HAS message "
                              f"IOD {iode} and observation epoch {obs_epoch.isot}. Use 'gnss_clean_orbit' remover.") 
               
                nearest_idx = self.get_nearest_idx(idx, obs_epoch.gps.mjd, time_key, positive, sat)
                brdc_idx.append(idx.nonzero()[0][nearest_idx])
            
        else:    

            # Determine broadcast ephemeris block index for a given satellite and observation epoch
            for sat, obs_epoch in zip(dset.satellite, dset[time]): # TODO: Is _add_time needed for dset[time]?
 
                idx = self.dset_edit.filter(satellite=sat)
               
                if np.any(idx) == False:
                    log.fatal(f"No valid broadcast navigation message could be found for satellite {sat} and "
                              f"observation epoch {obs_epoch.isot}. Use 'gnss_clean_orbit' remover.")   
                nearest_idx = self.get_nearest_idx(idx, obs_epoch.gps.mjd, time_key, positive, sat)
                brdc_idx.append(idx.nonzero()[0][nearest_idx])

        return brdc_idx
    

    def _get_corrected_broadcast_ephemeris(
        self, t_sat_gpsweek: float, t_sat_gpssec: float, idx: int, sys: str
    ) -> Dict[str, float]:
        """Apply correction for broadcast ephemeris for given time tk

        Following equations are based on Table 20-IV. in :cite:`is-gps-200h`.

        Args:
            t_sat_gpsweek (float):   GPS week of satellite transmission time.
            t_sat_gpssec (float):    GPS seconds of satellite transmission.
            idx (int):               Index for broadcast ephemeris dataset valid for observation time of receiver
            sys (str):               GNSS identifier

        Returns:
            dict:   Selected and prepared broadcast ephemeris dictionary with following entries:

        | Keys           | Unit  | Description
        | :------------- | :---- | :-------------------------------------------------------- |
        | a              | m     | Semimajor axis                                            |
        | E              | rad   | Eccentric anomaly                                         |
        | i              | rad   | Inclination                                               |
        | lambda_        | rad   | Instantaneous Greenwich longitude of the ascending node   |
        | n              | rad/s | Corrected mean motion                                     |
        | r              | m     | Orbit radius                                              |
        | tk             | s     | Eclapsed time referred to ephemeris reference epoch       |
        | u              | rad   | Argument of latitude                                      |
        | vega           | rad   | True anomaly                                              |

        """
        # BUG: Use of GPSSEC does not work for GPS WEEK crossovers. MJD * Unit.day2second() would a better solution. The
        #     problem is that use of GPSSEC compared to MJD * Unit.day2second() is not consistent!!!!
        gpsweek_diff = (t_sat_gpsweek - self._add_dim(self.dset_edit.toe.gps.gps_ws.week)[idx]) * Unit.week2second
        toe = self._add_dim(self.dset_edit.toe.gps.gps_ws.seconds)[idx] # Ephemeris reference epoch in [s]
        tk = t_sat_gpssec - toe + gpsweek_diff  # Eclapsed time referred to ephemeris reference epoch in [s]

        # Determine corrected Keplerian elements
        e = self.dset_edit.e[idx]  # Eccentricity of the orbit
        a = self.dset_edit.sqrt_a[idx] ** 2  # Semimajor axis in [m]
        omega = self.dset_edit.omega[idx]  # Argument of perigee in [rad]
        n0 = np.sqrt(GM[sys] / a ** 3)  # Mean motion of Keplerian orbit in [rad/s]
        n = n0 + self.dset_edit.delta_n[idx]  # Corrected mean motion in [rad/s]

        M = self.dset_edit.m0[idx] + n * tk  # Mean anomaly in [rad]
        E = self._get_eccentric_anomaly(M, e)  # Eccentric anomaly in [rad]
        vega = np.arctan2(np.sqrt(1.0 - e ** 2) * np.sin(E), np.cos(E) - e)  # True anomaly in [rad]
        u0 = omega + vega  # Initial argument of latitude
        lambda_ = self.dset_edit.Omega[idx] + (self.dset_edit.Omega_dot[idx] - OMEGA[sys]) * tk - OMEGA[sys] * toe
        # Instantaneous Greenwich longitude of the ascending node in [rad]

        # Determine argument of latitude, orbit radius and inclination
        sin2u0 = np.sin(2 * u0)
        cos2u0 = np.cos(2 * u0)
        du = self.dset_edit.cus[idx] * sin2u0 + self.dset_edit.cuc[idx] * cos2u0
        dr = self.dset_edit.crs[idx] * sin2u0 + self.dset_edit.crc[idx] * cos2u0
        di = self.dset_edit.cis[idx] * sin2u0 + self.dset_edit.cic[idx] * cos2u0
        u = u0 + du  # Argument of latitude in [rad]
        r = a * (1.0 - e * np.cos(E)) + dr  # Orbit radius in [m]
        i = self.dset_edit.i0[idx] + self.dset_edit.idot[idx] * tk + di  # Inclination in [rad]

        return {"a": a, "E": E, "i": i, "lambda_": lambda_, "n": n, "r": r, "tk": tk, "u": u, "vega": vega}


    @staticmethod
    def _get_eccentric_anomaly(M, e):
        r""" Computes the eccentric anomaly for elliptic orbits

        Newton's method used for determination of eccentric anomaly E as shown in Eq. 2.42 in :cite:`montenbruck2012`.

        Args:
            M (float):        Mean anomaly in [rad].
            e (float):        Eccentricity of the orbit in the interval [0, 1].

        Returns:
            Float:        Eccentric anomaly in [rad].

        Example:
            >>> _get_eccentric_anomaly(0.00817235075205, 0.00223578442819)
            0.0081906631043202408
        """
        num_iter = 0
        f = 1
        max_iter = 10
        limit = 10e-14  # TODO: Should limit be related to machine accuracy np.finfo(float)?

        # For small eccentriciy (e < 0.8) is it enough to start with E = M. For high eccentric orbits (e > 0.8) the
        # iteration should start with E = PI to avoid convergence problems during the iteration (see p. 24 in [1]).
        if e < 0.8:
            E = M
        else:
            E = np.pi

        while np.fabs(f) > limit:
            f = E - e * np.sin(E) - M
            E = E - f / (1.0 - e * np.cos(E))
            num_iter += 1

            if num_iter == max_iter:
                log.fatal(f"Convergence problem by determination of eccentric anomaly (max_iter = {max_iter})")

        return E
    

    def _get_relativistic_clock_correction(self, t_sat_gpsweek, t_sat_gpssec, idx, sys):
        """Determine relativistic clock correction due to orbit eccentricity for given observation epoch (satellite
        transmission time) based on Section 20.3.3.3.3.1 in :cite:`is-gps-200h`.

        Args:
            t_sat_gpsweek (float):  GPS week of satellite transmission time.
            t_sat_gpssec (float):   GPS seconds of satellite transmission.
            idx (int):              Index for broadcast ephemeris dataset valid for observation time of receiver
            sys (str):              GNSS identifier

        Returns:
            float64:      Relativistic orbit eccentricity correction in [m]
        """
        # Correct broadcast ephemeris
        bdict = self._get_corrected_broadcast_ephemeris(t_sat_gpsweek, t_sat_gpssec, idx, sys)

        # Compute relativistic orbit eccentricity in [m]
        return -2 / constant.c * np.sqrt(bdict["a"] * GM[sys]) * self.dset_edit.e[idx] * np.sin(bdict["E"])   


    def _get_satellite_position_vector(self, bdict):
        """Determine satellite position vector in Earth centered Earth fixed (ECEF) coordinate system.

        Following equations are based on Table 20-IV. in :cite:`is-gps-200h`.

        Args:
            bdict (dict):       Selected and prepared broadcast ephemeris dictionary with following entries

        | Keys           | Unit  | Description                                               |
        | :------------- | :---- | :-------------------------------------------------------- |
        | a              | m     | Semimajor axis                                            |
        | e              |       | Eccentricity of the orbit                                 |
        | E              | rad   | Eccentric anomaly                                         |
        | i              | rad   | Inclination                                               |
        | lambda_        | rad   | Instantaneous Greenwich longitude of the ascending node   |
        | n              | rad/s | Corrected mean motion                                     |
        | r              | m     | Orbit radius                                              |
        | tk             | s     | Eclapsed time referred to ephemeris reference epoch       |
        | u              | rad   | Argument of latitude                                      |
        | vega           | rad   | True anomaly                                              |

        Returns:
            dict:    Following entries are added to broadcast ephemeris dictionary

        | Keys           | Unit  | Description                                                                        |
        | :------------- | :---- | :--------------------------------------------------------------------------------- |
        | r_ecef         | m     | 3-dimensional numpy array with satellite position in Earth centered Earth fixed    |
        |                |       | (ECEF) coordinate system                                                           |
        | r_orb          | m     | 3-dimensional numpy array with satellite position in orbital coordinate system     |
        """
        # TODO: What should be done for GEO satellites (e.g. see in gLAB function getPositionBRDC() in model.c)?

        # Transformation from spherical to cartesian orbital coordinate system
        r_orb = np.array([bdict["r"] * np.cos(bdict["u"]), bdict["r"] * np.sin(bdict["u"]), 0.0])

        # Transformation from cartesian orbital to Earth centered Earth fixed (ECEF) geocentric equatorial coordinate
        # system
        r_ecef = rotation.R3(-bdict["lambda_"]).dot(rotation.R1(-bdict["i"])).dot(r_orb.T)

        bdict.update({"r_ecef": r_ecef, "r_orb": r_orb})


    def _get_satellite_position_velocity(self, t_sat_gpsweek, t_sat_gpssec, idx, sys):
        """Determine satellite position and velocity vector based on broadcast ephemeris

        Position and velocity is determined for given observation epoch (satellite transmission time) in the Earth-
        centered, Earth-fixed (ECEF) coordinate system.

        Args:
            t_sat_gpsweek (float):  GPS week of satellite transmission time.
            t_sat_gpssec (float):   GPS seconds of satellite transmission.
            idx (int):              Index for broadcast ephemeris dataset valid for observation time of receiver.
            sys (str):              GNSS identifier

        Returns:
            tuple:         with following elements

        | Elements       | Description                                                 |
        | :------------- | :---------------------------------------------------------- |
        | r_ecef         | Satellite position vector in ECEF coordinate system in [m]  |
        | v_ecef         | Satellite velocity vector in ECEF coordinate system in [m]  |
        """
        # TODO: get_row() from dataset based on idx

        # Correct broadcast ephemeris
        bdict = self._get_corrected_broadcast_ephemeris(t_sat_gpsweek, t_sat_gpssec, idx, sys)

        # Compute satellite position vector
        self._get_satellite_position_vector(bdict)

        # Compute satellite velocity vector
        self._get_satellite_velocity_vector(idx, bdict, sys)

        return bdict["r_ecef"], bdict["v_ecef"]


    def _get_satellite_velocity_vector(self, idx, bdict, sys):
        """Determine satellite velocity vector in Earth centered Earth fixed (ECEF) coordinate system.

        Following equations are described in :cite:`remondi2004`.

        Args:
            idx (int):          Index for broadcast ephemeris dataset valid for observation time of receiver.
            bdict (dict):       Selected and prepared broadcast ephemeris dictionary with following entries


        | Keys           | Unit  | Description
        | :------------- | :---- | :--------------------------------------------------------------------------------- |
        | a              | m     | Semimajor axis                                                                     |
        | E              | rad   | Eccentric anomaly                                                                  |
        | i              | rad   | Inclination                                                                        |
        | lambda_        | rad   | Instantaneous Greenwich longitude of the ascending node                            |
        | n              | rad/s | Corrected mean motion                                                              |
        | r              | m     | Orbit radius                                                                       |
        | r_orb          | m     | 3-dimensional numpy array with satellite position vector in orbital coordinate     |
        |                |       | system                                                                             |
        | r_ecef         | m     | 3-dimensional numpy array with satellite position vector in Earth centered         |
        |                |       | Earth-fixed coordinate system                                                      |
        | tk             | s     | Eclapsed time referred to ephemeris reference epoch                                |
        | u              | rad   | Argument of latitude                                                               |
        | vega           | rad   | True anomaly                                                                       |

        Returns:
            dict:   Following entries are added to broadcast ephemeris dictionary

        | Keys           | Unit  | Description                                                                        |
        | :------------- | :---- | :--------------------------------------------------------------------------------- |
        | v_ecef         | m     | 3-dimensional numpy array with satellite velocity vector in Earth centered         |
        |                |       | Earth-fixed coordinate system                                                      |
        | v_orb          | m     | 3-dimensional numpy array with satellite velocity vector in orbital coordinate     |
        |                |       | system                                                                             |
        """
        # Determine time derivatives of Keplerian elements
        M_dot = bdict["n"]  # Time derivative of mean anomaly in [rad/s]
        E_dot = M_dot / (
            1.0 - self.dset_edit.e[idx] * np.cos(bdict["E"])
        )  # Time derivative of eccentric anomaly in [rad/s]
        v_dot = (
            (1.0 + self.dset_edit.e[idx] * np.cos(bdict["vega"]))
            * E_dot
            * np.sin(bdict["E"])
            / ((1.0 - self.dset_edit.e[idx] * np.cos(bdict["E"])) * np.sin(bdict["vega"]))
        )
        # Time derivative of true anomaly in [rad/s]
        u0_dot = v_dot  # Time derivative of initial argument of latitude in [rad/s]
        lambda_dot = self.dset_edit.Omega_dot[idx] - OMEGA[sys]  # Time derivative of instantaneous Greenwich longitude
        # of the ascending node in [rad]

        # Determine time derivatives of argument of latitude, orbit radius and inclination
        sin2u = np.sin(2 * bdict["u"])
        cos2u = np.cos(2 * bdict["u"])
        du_dot = 2 * u0_dot * (self.dset_edit.cus[idx] * cos2u - self.dset_edit.cuc[idx] * sin2u)
        dr_dot = 2 * u0_dot * (self.dset_edit.crs[idx] * cos2u - self.dset_edit.crc[idx] * sin2u)
        di_dot = 2 * u0_dot * (self.dset_edit.cis[idx] * cos2u - self.dset_edit.cic[idx] * sin2u)
        u_dot = u0_dot + du_dot  # Time derivative of argument of latitude in [rad/s]
        r_dot = bdict["a"] * self.dset_edit.e[idx] * E_dot * np.sin(bdict["E"]) + dr_dot
        # Time derivative of orbit radius in [m/s]
        i_dot = self.dset_edit.idot[idx] + di_dot  # Time derivative of inclination in [rad/s]

        # Determine satellite velocity vector in the orbital plane coordinate system
        v_orb = np.array(
            [
                r_dot * np.cos(bdict["u"]) - bdict["r"] * u_dot * np.sin(bdict["u"]),
                r_dot * np.sin(bdict["u"]) + bdict["r"] * u_dot * np.cos(bdict["u"]),
                0.0,
            ]
        )

        # Satellite velocity vector in Earth centered Earth-fixed coordinate system
        x_dot = (
            v_orb[0] * np.cos(bdict["lambda_"])
            - v_orb[1] * np.cos(bdict["i"]) * np.sin(bdict["lambda_"])
            + bdict["r_orb"][1] * i_dot * np.sin(bdict["i"]) * np.sin(bdict["lambda_"])
            - bdict["r_ecef"][1] * lambda_dot
        )

        y_dot = (
            v_orb[0] * np.sin(bdict["lambda_"])
            + v_orb[1] * np.cos(bdict["i"]) * np.cos(bdict["lambda_"])
            - bdict["r_orb"][1] * i_dot * np.sin(bdict["i"]) * np.cos(bdict["lambda_"])
            + bdict["r_ecef"][0] * lambda_dot
        )

        z_dot = v_orb[1] * np.sin(bdict["i"]) + bdict["r_orb"][1] * i_dot * np.cos(bdict["i"])

        v_ecef = np.array([x_dot, y_dot, z_dot])

        bdict.update({"v_ecef": v_ecef, "v_orb": v_orb})


