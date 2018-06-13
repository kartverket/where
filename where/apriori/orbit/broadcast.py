"""Handling of GNSS broadcast orbits

Description:
------------
The module includes a class for handling apriori GNSS broadcast orbits.

Example:
    from where import apriori

    # Get broadcast ephemeris object
    brdc = apriori.get('orbit', rundate=rundate, time=time, satellite=satellite, system=system, station=station,
                       apriori_orbit='broadcast')

    # Write calculated Dataset to file
    brdc.dset.write()



"""

# Standard library imports
import re

# External library imports
import numpy as np
import pandas as pd

# Where imports
from where import apriori
from where import cleaners
from where import data
from where import parsers
from where.apriori import orbit
from where.lib import config
from where.lib import constant
from where.lib import files
from where.lib import gnss
from where.lib import log
from where.lib import mathp
from where.lib import plugins
from where.lib import rotation
from where.lib.unit import unit

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
        satellite (tuple):      Satellite number together with GNSS system identifier given for each observation epoch
        system (tuple):         GNSS system identifier given for each observation epoch
        time (Time):            Observation epochs

    Methods:
        relativistic_clock_correction():  Determine relativistic clock correction due to orbit eccentricity
        satellite_clock_correction():  Determine satellite clock correction
        unhealthy_satellites(): Get unhealthy satellites based on RINEX navigation file information
        update_time():          Update which time epochs to calculate orbit for

        _calculate():           Calculate broadcast ephemeris and satellite clock correction for given observation
                                epochs
        _edit():                Edit RINEX navigation file data and save it in a Dataset
        _read():                Read RINEX navigation file data and save it in a Dataset
    """
    name = "broadcast"

    def __init__(self, rundate, time=None, satellite=None, system=None, station=None, file_key=None):
        """Set up a new BroadcastOrbit object, does not parse any data

        TODO: Remove dependency on rundate, use time to read correct files. (What to do with dataset?)

        Args:
            rundate (date):     Date of model run.
            time (Time):        Time epochs at the satellite for which to calculate the apriori orbit.
            satellite (list):   Strings with names of satellites.
            system (list):      Strings with codes of systems (G, E, R, etc.).
            station (str):      4 character station identifier.
            file_key (str):     Key to the broadcast orbit file defined in files.conf file.
        """
        super().__init__(rundate=rundate, time=time, satellite=satellite, system=system)
        self.file_key = "gnss_rinex_nav_{system}" if file_key is None else file_key

        # TODO hjegei: Should it be not enough to 'station' in _dset_raw?
        self._dset_raw.vars["station"] = station.lower()
        self._dset_raw.vars["STATION"] = self._dset_raw.vars["station"].upper()
        self._dset_edit.vars["station"] = station.lower()
        self._dset_edit.vars["STATION"] = self._dset_raw.vars["station"].upper()

    def _read(self, dset_raw):
        """Read RINEX navigation file data and save it in a Dataset

        One RINEX navigation file is normally written for each GNSS. The navigation file extension depends on the GNSS
        (GPS: *.n, Galileo: *.l, GLONASS: *.g, ...). Therefore we have defined for each GNSS the navigation file
        name in the Where configuration file `files.conf`. In addition mixed navigation files exists, which includes
        navigation messages of different GNSS. We use following file keys in `files.conf`:
            =========  ==================
             System     File key
            =========  ==================
             Galileo    gnss_rinex_nav_E
             GLONASS    gnss_rinex_nav_R
             GPS        gnss_rinex_nav_G
             Mixed      gnss_rinex_nav_M
            =========  ==================

        Depending on the configuration options `systems` and `use_mixed_brdc_file` following navigation files are read:

         ======================  ==================  =======================================
          Option                  File key            What kind of navigation file is read?
         ======================  ==================  =======================================
          systems = G             gnss_rinex_nav_G    Only the GPS navigation file
          systems = G E           gnss_rinex_nav_G    GPS and Galileo navigation files
                                  gnss_rinex_nav_E
          use_mixed_brdc_file     gnss_rinex_nav_M    Mixed GNSS navigation file
         ======================  ==================  =======================================

        Args:
            dset_raw (Dataset):     Dataset representing raw data from RINEX navigation file
        """
        use_mixed_brdc_file = config.tech.get("use_mixed_brdc_file", default=False).bool
        systems = {"M"} if use_mixed_brdc_file == True else set(self.system)
        file_paths = list()

        for sys in systems:

            # Generate temporary Dataset with orbit file data
            dset_temp = data.Dataset(dset_raw.rundate, dset_raw.vars["tech"], "temporary", "", 0, empty=True)
            parser = parsers.parse(
                file_key=self.file_key.format(system=sys),
                rundate=dset_raw.rundate,
                file_vars=dict(dset_raw.vars, file_key=self.file_key.format(system=sys)),
            )
            parser.write_to_dataset(dset_temp)
            file_paths.append(str(parser.file_path))

            # Extend Dataset dset_raw with temporary Dataset
            if dset_raw.num_obs:

                # Merge meta data information
                # TODO: Handle meta data information correctly. Meta data information based on different GNSS navigation
                #       message files has to be merged correctly together. What to do with 'sat_sys' and 'leap_seconds'?
                for date in dset_temp.meta.keys():
                    for key in ["iono_para", "time_sys_corr"]:
                        dset_temp.meta[date][key].update(dset_raw.meta[date][key])

                dset_raw.extend(dset_temp)
            else:
                dset_raw.copy_from(dset_temp)
            dset_raw.add_to_meta("parser", "file_path", file_paths)

        return dset_raw

    def _edit(self, dset_edit):
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
        nav = self.dset_raw.as_dataframe()

        # Filter frame
        nav_filtered = nav.sort_values(by=["satellite", "time", "transmission_time"])

        # Remove duplicated navigation message entries (IODEs)
        # TODO: transmission_time should also be used in filtering!!! Maybe it should be a configuration option, how
        # to filter duplicated epochs. Keep first and last.
        idx = nav_filtered.duplicated(subset=["satellite", "time", "iode", "nav_type"], keep="first")
        with pd.option_context("display.max_rows", None, "display.max_columns", 5):
            log.debug(
                "Remove following duplicated navigation message entries: \n{}",
                nav_filtered[["satellite", "time", "iode", "nav_type"]][idx],
            )
        nav_filtered = nav_filtered.drop_duplicates(subset=["satellite", "time", "iode", "nav_type"], keep="first")

        # Remove navigation message types, which are not needed
        nav_type = config.tech.get("navigation_message_type", default="").dict
        if nav_type:
            for sys, type_ in nav_type.items():
                sys = sys[0] if len(sys) > 0 else sys
                if type_ == "INAV":
                    type_ = ["INAV_E1", "INAV_E5b", "INAV_E1E5b"]
                elif type_ == "FNAV":
                    type_ = ["FNAV_E5a"]

                remove_nav_type = nav_filtered.query("system == @sys and nav_type != @type_")
                if not remove_nav_type.empty:
                    log.info(
                        "Remove '{}' navigation messages of GNSS '{}'.",
                        ", ".join(list(set(remove_nav_type.nav_type))),
                        sys,
                    )
                    nav_filtered = pd.concat([nav_filtered, remove_nav_type]).drop_duplicates(keep=False)

        # TODO hjegei: Possible future ...
        # dset_edit.copy_from(self.dset_raw)
        # dset_edit.reorder(nav_filtered.index.values)

        # Convert edited fields to Dataset
        nav_np = nav_filtered.as_matrix()
        fields = nav_filtered.columns

        dset_edit.vars["orbit"] = self.name
        dset_edit.num_obs = nav_filtered.shape[0]
        dset_edit.meta.update(self.dset_raw.meta)

        for idx, field in enumerate(fields):
            if field in ["time", "transmission_time", "toe"]:
                dset_edit.add_time(field, val=nav_np[:, idx], scale="gps", format="datetime")
            elif field in ["nav_type", "satellite", "system"]:
                dset_edit.add_text(field, val=nav_np[:, idx])
            else:
                dset_edit.add_float(field, val=nav_np[:, idx].astype(float))

        return dset_edit

    def _calculate(self, dset):
        """Calculate broadcast ephemeris and satellite clock correction for given observation epochs

        Args:
            dset (Dataset): Dataset representing calculated broadcast ephemeris with following fields:

        ========================  ===============  =======  ========================================================
         Field                     Type             Unit     Description
        ========================  ===============  =======  ========================================================
         gnss_satellite_clock     numpy.ndarray     m       Satellite clock correction
         gnss_relativistic_clock  numpy.ndarray     m       Relativistic clock correction due to orbit eccentricity
         sat_posvel               PosVelTable       m       Satellite position and velocity
         satellite                numpy.ndarray             Satellite numbers
         system                   numpy.ndarray             GNSS identifiers
         time                     TimeTable                 Observation epochs
         used_iode                numpy.ndarray             IODE of selected broadcast ephemeris block
         used_transmission_time   TimeTable                 Transmission time of selected broadcast ephemeris block
         used_toe                 TimeTable                 Time of ephemeris (TOE) of selected broadcast ephemeris
                                                            block
        =======================  ===============  =======  ========================================================
        """
        not_available_sat = set(self.satellite) - set(self.dset_edit.satellite)
        if not_available_sat:
            log.fatal(
                "Satellites {} are not given in broadcast ephemeris file {}. Add satellite to 'ignore_satellite'"
                "option in configuration file.",
                ", ".join(sorted(not_available_sat)),
                ", ".join(self.dset_edit.meta["parser"]["file_path"]),
            )

        not_implemented_sys = set(self.system) - set("EG")
        if not_implemented_sys:
            log.fatal(
                "At the moment Where can provide broadcast ephemeris for GNSS 'E' and 'G', but not " "for {}.",
                not_implemented_sys,
            )

        log.info(
            "Calculating satellite position/velocity (broadcast) based on RINEX navigation file {}.",
            ", ".join(self.dset_edit.meta["parser"]["file_path"]),
        )

        dset.vars["orbit"] = self.name
        dset.num_obs = len(self.time)
        dset.add_time("time", val=self.time, scale=self.time.scale)
        dset.add_text("satellite", val=self.satellite)
        dset.add_text("system", val=self.system)

        sat_pos = np.zeros((dset.num_obs, 3))
        sat_vel = np.zeros((dset.num_obs, 3))

        # Get correct navigation block for given observations times by determining the indices to broadcast
        # ephemeris Dataset
        dset_brdc_idx = self._get_brdc_block_idx()

        # Loop over all observations
        # TODO: Generation of vectorized solution, if possible?

        # BUG: Use of GPSSEC does not work for GPS WEEK crossovers. MJD * unit.day2second() would a better solution. The
        #     problem is that use of GPSSEC compared to MJD * unit.day2second() is not consistent!!!!
        for obs_idx, (time_gpsweek, time_gpssec, brdc_idx, sys) in enumerate(
            zip(self.time.gps.gpsweek, self.time.gps.gpssec, dset_brdc_idx, self.system)
        ):

            # TODO: get_row() function needed for brdc -> brdc.get_row(kk)
            sat_pos[obs_idx], sat_vel[obs_idx] = self._get_satellite_position_velocity(
                time_gpsweek, time_gpssec, brdc_idx, sys
            )

            # +DEBUG
            # print("DEBUG: {} obs_idx: {:>5d} brdc_idx: {:>5d} toc: {:>5.0f} {:>6.0f} toe: {:>6.0f} trans_time: {:>6.0f}"
            #      " tk: {:>16.10f} iode: {:>3d} sqrt_a: {:>17.10f} sat_pos: {:>21.10f} {:>21.10f} {:>21.10f} "
            #      "sat_vel: {:>17.10f} {:>17.10f} {:>17.10f} sat_clk_bias: {:>17.10f}, sat_clk_drft: {:>17.10f} "
            #      ''.format(self.dset_edit.satellite[brdc_idx], obs_idx, brdc_idx,
            #                dset.time.gps.jd_frac[obs_idx] * 86400,
            #                dset.time.gps.gpssec[obs_idx],
            #                self.dset_edit.toe.gps.gpssec[brdc_idx],
            #                self.dset_edit.transmission_time.gps.gpssec[brdc_idx],
            #                dset.time.gps.jd_frac[obs_idx]-self.dset_edit.toe.gps.gpssec[brdc_idx],
            #                int(self.dset_edit.iode[brdc_idx]),
            #                self.dset_edit.sqrt_a[brdc_idx],
            #                sat_pos[obs_idx][0], sat_pos[obs_idx][1], sat_pos[obs_idx][2],
            #                sat_vel[obs_idx][0], sat_vel[obs_idx][1], sat_vel[obs_idx][2],
            #                self.dset_edit.sat_clock_bias[brdc_idx],
            #                self.dset_edit.sat_clock_drift[brdc_idx],))
            # -DEBUG

        # Add information about choosed broadcast ephemeris block
        dset.add_float("used_iode", val=self.dset_edit.iode[dset_brdc_idx])
        dset.add_time(
            "used_transmission_time",
            val=self.dset_edit.transmission_time[dset_brdc_idx],
            scale=self.dset_edit.transmission_time.scale,
        )
        dset.add_time("used_toe", val=self.dset_edit.toe[dset_brdc_idx], scale=self.dset_edit.toe.scale)

        # Add information about group delays
        for field in ["bgd_e1_e5a", "bgd_e1_e5b", "tgd", "tgd_b1_b3", "tgd_b2_b3"]:
            if field in self.dset_edit.fields:
                dset.add_float(field, val=self.dset_edit[field][dset_brdc_idx])

        # Add satellite clock correction to Dataset
        dset.add_float("gnss_satellite_clock", val=self.satellite_clock_correction(), unit="meter")

        # Add satellite position and velocity to Dataset
        dset.add_posvel("sat_posvel", time="time", itrs=np.hstack((sat_pos, sat_vel)))

        # Add relativistic clock correction to Dataset
        dset.add_float("gnss_relativistic_clock", val=self.relativistic_clock_correction(), unit="meter")

        # +DEBUG
        # for num_obs  in range(0, dset.num_obs):
        #    print('DEBUG: ', dset.satellite[num_obs],
        #                     dset.time.gps.datetime[num_obs],
        #                     dset.time.gps.mjd_frac[num_obs]*24*3600,
        #                     dset.gnss_satellite_clock[num_obs],
        #                     dset.gnss_relativistic_clock[num_obs])
        # -DEBUG

        return dset

    def satellite_clock_correction(self):
        """Determine satellite clock correction

        The satellite clock correction is based on Section 20.3.3.3.3.1 in :cite:`is-gps-200h`.

        Returns:
            numpy.ndarray:    GNSS satellite clock corrections for each observation in [m] (Note: without relativistic
                              orbit eccentricity correction)
        """
        # Get correct navigation block for given observations times by determining the indices to broadcast ephemeris
        # Dataset
        dset_brdc_idx = self._get_brdc_block_idx()

        # Elapsed time referred to clock data reference epoch toc in [s]
        # BUG: Use of GPSSEC does not work for GPS WEEK crossovers. MJD * unit.day2second() would a better solution. The
        #     problem is that use of GPSSEC compared to MJD * unit.day2second() is not consistent!!!!
        gpsweek_diff = (self.dset.time.gps.gpsweek - self.dset_edit.time.gps.gpsweek[dset_brdc_idx]) * unit.week2second
        tk = self.dset.time.gps.gpssec - self.dset_edit.time.gps.gpssec[dset_brdc_idx] + gpsweek_diff

        return (
            self.dset_edit.sat_clock_bias[dset_brdc_idx]
            + self.dset_edit.sat_clock_drift[dset_brdc_idx]
            * tk
            + self.dset_edit.sat_clock_drift_rate[dset_brdc_idx]
            * tk
            ** 2
        ) * constant.c

    def unhealthy_satellites(self):
        """Get unhealthy satellites based on RINEX navigation file information

        A satellite is set to unhealthy, if in one of the given broadcast ephemeris blocks the condition
        'self.dset_edit.sv_health > 0' is valid.

        Returns:
            list:    Unhealthy satellites
        """
        unhealthy_satellites = list()
        for satellite in self.dset_edit.unique("satellite"):
            idx = self.dset_edit.filter(satellite=satellite)
            if np.all(self.dset_edit.sv_health[idx] > 0):
                unhealthy_satellites.append(satellite)

        return unhealthy_satellites

    def _get_brdc_block_idx(self):
        """Get GNSS broadcast ephemeris block indices for given observation epochs

        The indices relate the observation epoch to the correct set of broadcast ephemeris. First the time difference
        between the observation epoch and a selected time is calculated to determine the correct broadcast ephemeris
        block. The seleted time can be either the navigation epoch (time of clock (TOC)), the time of ephemeris (TOE)
        or the transmission time. Afterwards the broastcast block is selected with the smallest time difference.
        Following option can be choosen for configuration file option 'brdc_block_nearest_to':

        ==============================  ================================================================================
         Option                           Description
        ==============================  ================================================================================
          toc                            Broadcast block for given observation epoch is selected nearest to navigation
                                         epoch (time of clock (TOC)).
          toc:positive                   Same as 'toc' option, but the difference between observation epoch and TOC has
                                         to be positive.
          toe                            Broadcast block for given observation epoch is selected nearest to time of
                                         ephemeris (TOE).
          toe:positive                   Same as 'toe' option, but the difference between observation epoch and TOE has
                                         to be positive.
          transmission_time              Broadcast block for given observation epoch is selected nearest to transmission
                                         time.
          transmission_time:positive     Same as 'transmission_time' option, but the difference between observation
                                         epoch and transmission time has to be positive.
        =============================  =================================================================================

        Returns:
            List: Broadcast ephemeris block indices for given observation epochs.
        """

        brdc_block_nearest_to_options = [
            "toc", "toc:positive", "toe", "toe:positive", "transmission_time", "transmission_time:positive"
        ]
        brdc_idx = list()

        # Get configuration option
        brdc_block_nearest_to = config.tech.get("brdc_block_nearest_to", default="toe:positive").str.rsplit(":", 1)
        if ":".join(brdc_block_nearest_to) not in brdc_block_nearest_to_options:
            log.fatal(
                "Unknown value '{}' for configuration file option 'brdc_block_nearest_to'. Following values can "
                "be selected: {}.",
                ":".join(brdc_block_nearest_to),
                ", ".join(brdc_block_nearest_to_options),
            )

        time_key = brdc_block_nearest_to[0]
        positive = True if "positive" in brdc_block_nearest_to else False

        log.debug("Broadcast block is selected nearest to '{}{}' time.", "+" if positive else "+/-", time_key)

        # Check if broadcast orbits are available
        not_available_sat = set(self.satellite) - set(self.dset_edit.satellite)
        if not_available_sat:
            log.fatal(
                "Following satellites are not given in apriori broadcast orbit file {}: {}",
                ", ".join(self.dset_edit.meta["parser"]["file_path"]),
                ", ".join(sorted(not_available_sat)),
            )

        # Determine broadcast ephemeris block index for a given satellite and observation epoch
        for sat, time in zip(self.satellite, self.time):

            idx = self.dset_edit.filter(satellite=sat)
            diff = time.gps.mjd - self.dset_edit[time_key].gps.mjd[idx]
            if positive:
                nearest_idx = np.array([99999 if v < 0 else v for v in diff]).argmin()
            else:
                nearest_idx = np.array([abs(diff)]).argmin()

            brdc_idx.append(idx.nonzero()[0][nearest_idx])

        return brdc_idx

    def _get_corrected_broadcast_ephemeris(self, t_sat_gpsweek, t_sat_gpssec, idx, sys):
        """Apply correction for broadcast ephemeris for given time tk

        Following equations are based on Table 20-IV. in :cite:`is-gps-200h`.

        Args:
            t_sat_gpsweek (float):   GPS week of satellite transmission time.
            t_sat_gpssec (float):    GPS seconds of satellite transmission.
            idx (int):               Index for broadcast ephemeris dataset valid for observation time of receiver
            sys (str):               GNSS identifier

        Returns:
            dict:   Selected and prepared broadcast ephemeris dictionary with following entries:

        ===============  ======  =============================================================
         Keys             Unit    Description
        ===============  ======  =============================================================
         a                m       Semimajor axis
         E                rad     Eccentric anomaly
         i                rad     Inclination
         lambda_          rad     Instantaneous Greenwich longitude of the ascending node
         n                rad/s   Corrected mean motion
         r                m       Orbit radius
         tk               s       Eclapsed time referred to ephemeris reference epoch
         u                rad     Argument of latitude
         vega             rad     True anomaly
        ===============  ======  =============================================================
        """
        toe = self.dset_edit.toe.gps.gpssec[idx]  # Ephemeris reference epoch in [s]
        # BUG: Use of GPSSEC does not work for GPS WEEK crossovers. MJD * unit.day2second() would a better solution. The
        #     problem is that use of GPSSEC compared to MJD * unit.day2second() is not consistent!!!!
        gpsweek_diff = (t_sat_gpsweek - self.dset_edit.toe.gps.gpsweek[idx]) * unit.week2second
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
                log.fatal("Convergence problem by determination of eccentric anomaly (max_iter = {}).", max_iter)

        return E

    def _get_satellite_position_vector(self, bdict):
        """Determine satellite position vector in Earth centered Earth fixed (ECEF) coordinate system.

        Following equations are based on Table 20-IV. in :cite:`is-gps-200h`.

        Args:
            bdict (dict):       Selected and prepared broadcast ephemeris dictionary with following entries

        ===============  ======  =============================================================
         Keys             Unit    Description
        ===============  ======  =============================================================
         a                m       Semimajor axis
         e                        Eccentricity of the orbit
         E                rad     Eccentric anomaly
         i                rad     Inclination
         lambda_          rad     Instantaneous Greenwich longitude of the ascending node
         n                rad/s   Corrected mean motion
         r                m       Orbit radius
         tk               s       Eclapsed time referred to ephemeris reference epoch
         u                rad     Argument of latitude
         vega             rad     True anomaly
        ===============  ======  =============================================================

        Returns:
            dict:    Following entries are added to broadcast ephemeris dictionary

        ===============  ======  ==================================================================================
         Keys             Unit    Description
        ===============  ======  ==================================================================================
         r_ecef           m       3-dimensional numpy array with satellite position in Earth centered Earth fixed (ECEF)
                                  coordinate system
         r_orb            m       3-dimensional numpy array with satellite position in orbital coordinate system
        ===============  ======  ==================================================================================
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

        ===============  ==================================================================================
         Elements         Description
        ===============  ==================================================================================
         r_ecef           Satellite position vector in ECEF coordinate system in [m]
         v_ecef           Satellite velocity vector in ECEF coordinate system in [m]
        ===============  ==================================================================================
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

        ===============  ======  ==================================================================================
         Keys             Unit    Description
        ===============  ======  ==================================================================================
         a                m       Semimajor axis
         E                rad     Eccentric anomaly
         i                rad     Inclination
         lambda_          rad     Instantaneous Greenwich longitude of the ascending node
         n                rad/s   Corrected mean motion
         r                m       Orbit radius
         r_orb            m       3-dimensional numpy array with satellite position vector in orbital coordinate
                                  system
         r_ecef          m        3-dimensional numpy array with satellite position vector in Earth centered Earth-fixed
                                  coordinate system
         tk               s       Eclapsed time referred to ephemeris reference epoch
         u                rad     Argument of latitude
         vega             rad     True anomaly
        ===============  ======  ==================================================================================

        Returns:
            dict:   Following entries are added to broadcast ephemeris dictionary

        ===============  ======  ==================================================================================
         Keys             Unit    Description
        ===============  ======  ==================================================================================
         v_ecef          m        3-dimensional numpy array with satellite velocity vector in Earth centered Earth-fixed
                                  coordinate system
         v_orb            m       3-dimensional numpy array with satellite velocity vector in orbital coordinate
                                  system
        ===============  ======  ==================================================================================
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
            v_orb[0]
            * np.cos(bdict["lambda_"])
            - v_orb[1]
            * np.cos(bdict["i"])
            * np.sin(bdict["lambda_"])
            + bdict["r_orb"][1]
            * i_dot
            * np.sin(bdict["i"])
            * np.sin(bdict["lambda_"])
            - bdict["r_ecef"][1]
            * lambda_dot
        )

        y_dot = (
            v_orb[0]
            * np.sin(bdict["lambda_"])
            + v_orb[1]
            * np.cos(bdict["i"])
            * np.cos(bdict["lambda_"])
            - bdict["r_orb"][1]
            * i_dot
            * np.sin(bdict["i"])
            * np.cos(bdict["lambda_"])
            + bdict["r_ecef"][0]
            * lambda_dot
        )

        z_dot = v_orb[1] * np.sin(bdict["i"]) + bdict["r_orb"][1] * i_dot * np.cos(bdict["i"])

        v_ecef = np.array([x_dot, y_dot, z_dot])

        bdict.update({"v_ecef": v_ecef, "v_orb": v_orb})

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
