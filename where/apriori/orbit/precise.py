"""Handling of GNSS precise orbits

Description:
------------
The module includes a class for handling apriori GNSS precise orbits.

Example:
    from where import apriori

    # Get broadcast ephemeris object
    precise = apriori.get('orbit', rundate=rundate, time=time, satellite=satellite, system=system, station=station,
                          apriori_orbit='precise')

    # Write calculated Dataset to file
    precise.dset.write()


$Revision: 15127 $
$Date: 2018-05-16 19:46:18 +0200 (Wed, 16 May 2018) $
$LastChangedBy: hjegei $

"""
# Standard library imports
from datetime import timedelta

# External library imports
from midgard.math import interpolation
import numpy as np
import pandas as pd
from scipy import interpolate

# Where imports
from where import apriori
from where import data
from where.apriori import orbit
from where.lib import config
from where.lib import files
from where.lib import log
from where.lib import mathp
from where.lib import plugins
from where import parsers


@plugins.register
class PreciseOrbit(orbit.AprioriOrbit):
    """A class for representing apriori precise orbits

    SP3 orbit files can be read. In addition GNSS satellite position, velocities, satellite clock correction and
    relativistic clock correction can be calculated for given observation epochs based on precise orbits.

    The satellite position is determined for each observation epoch by interpolating within the given SP3 orbit time
    entries. The satellite velocities are calculated based on satellite positions 0.5 second before and after the
    observation epoch.

    Attributes:
        day_offset (int):       Day offset used to calculate the number of days to read.
        dset (Dataset):         Dataset object, which includes satellite position and velocity, satellite clock
                                correction and relativistic clock correction for each observation epoch
        dset_raw (Dataset):     Dataset object, which includes precise orbits read from SP3 file
        file_key (str):         Key to the precise orbit file defined in files.conf file.
        file_path (pathlib.PosixPath):  File path to SP3 orbit file.
        name (str):             Apriori orbit name
        satellite (tuple):      Satellite number together with GNSS system identifier given for each observation epoch
        system (tuple):         GNSS system identifier given for each observation epoch
        time (Time):            Observation epochs

    Methods:
        satellite_clock_correction():  Determine satellite clock correction
        _calculate():           Calculate precise orbit and satellite clock correction for given observation epochs
        _edit():                Edit precise orbit file data and save it in a Dataset
        _read():                Read precise orbit file data and save it in a Dataset

    """
    name = "precise"

    def __init__(self, rundate, time, satellite, system=None, file_key=None, file_path=None, day_offset=1, **kwargs):
        """Set up a new PreciseOrbit object, does not parse any data

        TODO: Remove dependency on rundate, use time to read correct files. (What to do with dataset?)

        Args:
            rundate (date):     Date of model run.
            time (Time):        Time epochs at the satellite for which to calculate the apriori orbit.
            satellite (list):   Strings with names of satellites.
            system (list):      Strings with codes of systems (G, E, R, etc.).
            file_key (str):     Key to the precise orbit file that will be read.
            file_path (pathlib.PosixPath):  File path to SP3 orbit file.
            day_offset (int):   Day offset used to calculate the number of days to read.
        """
        super().__init__(rundate=rundate, time=time, satellite=satellite, system=system)
        self.file_key = "gnss_orbit_sp3" if file_key is None else file_key
        self.day_offset = day_offset
        self.file_path = file_path

    def _read(self, dset_raw):
        """Read SP3 orbit file data and save it in a Dataset

        In addition to the given date, we read data for the day before and after. This is needed to carry out correct
        orbit interpolation at the start and end of a day.

        TODO:
        How well fits the orbits from day to day? Is it necessary to align the orbits?

        Args:
            dset_raw (Dataset):   Dataset representing raw data from apriori orbit files
        """
        date_to_read = dset_raw.rundate - timedelta(days=self.day_offset)
        file_paths = list()

        # Loop over days to read
        while date_to_read <= dset_raw.rundate + timedelta(days=self.day_offset):
            if self.file_path is None:
                file_path = files.path(self.file_key, file_vars=config.date_vars(date_to_read))
            else:
                file_path = self.file_path

            log.debug("Parse precise orbit file {}.", file_path)

            # Generate temporary Dataset with orbit file data
            dset_temp = data.Dataset(
                rundate=date_to_read,
                tech=dset_raw.vars["tech"],
                stage="temporary",
                dataset_name="",
                dataset_id=0,
                empty=True,
            )
            parser = parsers.parse(parser_name="orbit_sp3", file_path=file_path, rundate=date_to_read)
            parser.write_to_dataset(dset_temp)
            file_paths.append(str(parser.file_path))

            # Extend Dataset dset_raw with temporary Dataset
            dset_raw.copy_from(dset_temp) if dset_raw.num_obs == 0 else dset_raw.extend(dset_temp)
            dset_raw.add_to_meta("parser", "file_path", file_paths)

            date_to_read += timedelta(days=1)

        return dset_raw

    def _edit(self, dset_edit):
        """Edit precise orbit data and save it in a Dataset

        First the precise orbits are sorted after the satellite and time. Afterwards duplicated precise orbits for a 
        satellite are removed, whereby the last accurance of the precise orbit entry is kept (Note: SP3 files include
        often at the end the first epoch of the next day, whereby satellite clock corrections are not given. Therefore
        we keep the last accurance to be sure to have satellite clock corrections values for the first epoch of a day). 

        Args:
            dset_edit (Dataset):     Dataset representing edited data from precise orbit file
        """

        # TODO: Would it make sense to move this functionality to a "remover"?

        # Generate Pandas frame for all precise orbit entries
        precise = self.dset_raw.as_dataframe()

        # Filter frame
        precise_filtered = precise.sort_values(by=["satellite", "time"])

        # Remove duplicated precise orbit entries
        idx = precise_filtered.duplicated(subset=["satellite", "time"], keep="last")
        with pd.option_context("display.max_rows", None, "display.max_columns", 5):
            log.debug(
                "Remove following duplicated precise orbit entries: \n{}", precise_filtered[["satellite", "time"]][idx]
            )
        precise_filtered = precise_filtered.drop_duplicates(subset=["satellite", "time"], keep="last")

        # TODO hjegei: Possible future ...
        # dset_edit.copy_from(self.dset_raw)
        # dset_edit.reorder(precise_filtered.index.values)

        # Convert edited fields to Dataset
        precise_np = precise_filtered.as_matrix()
        fields = precise_filtered.columns

        dset_edit.vars["orbit"] = self.name
        dset_edit.num_obs = precise_filtered.shape[0]
        dset_edit.meta.update(self.dset_raw.meta)

        for idx, field in enumerate(fields):
            if field in ["time"]:
                # Note: It is ensured time scale consistency by using time scale from self.dset_raw. This is important,
                #       because 'self.dset_raw.as_dataframe()' uses time scale, when self.dset_raw time object was
                #       initialized.
                dset_edit.add_time(field, val=precise_np[:, idx], scale=self.dset_raw.time.scale, format="datetime")
            elif field in ["satellite", "system"]:
                dset_edit.add_text(field, val=precise_np[:, idx])
            elif field in ["sat_clock_bias"]:
                dset_edit.add_float(field, val=precise_np[:, idx].astype(float))

        # Add satellite position values to Dataset
        sat_pos = np.stack(
            (
                precise_np[:, fields.get_loc("sat_pos_itrs_0")],
                precise_np[:, fields.get_loc("sat_pos_itrs_1")],
                precise_np[:, fields.get_loc("sat_pos_itrs_2")],
            ),
            axis=1,
        )
        dset_edit.add_position("sat_pos", time="time", itrs=sat_pos, unit="meter")

        return dset_edit

    def _calculate(self, dset):
        """Calculate precise orbits and satellite clock correction for given observation epochs

        The satellite position is determined for each observation epoch by interpolating within the given SP3 orbit time
        entries. The satellite velocities are calculated based on satellite positions 0.5 second before and after the
        observation epoch.

        Args:
            dset (Dataset): Dataset representing calculated precise orbits with following fields:

        ========================  ===============  =======  ========================================================
         Field                     Type             Unit     Description
        ========================  ===============  =======  ========================================================
         gnss_satellite_clock     numpy.ndarray     m       Satellite clock correction
         gnss_relativistic_clock  numpy.ndarray     m       Relativistic clock correction due to orbit eccentricity
         sat_posvel               PosVelTable       m       Satellite position and velocity
         satellite                numpy.ndarray             Satellite numbers
         system                   numpy.ndarray             GNSS identifiers
         time                     TimeTable                 Observation epochs
        =======================  ===============  =======  ========================================================
        """
        # Check if satellites are given in SP3 file
        # TODO: Another solution has to be found for satellites not given in SP3 file, e.g. use of broadcast
        #       ephemeris.
        not_available_sat = set(self.satellite) - set(self.dset_edit.satellite)
        if not_available_sat:
            log.fatal(
                "Satellites {} are not given in precise orbit file {}.",
                ", ".join(sorted(not_available_sat)),
                self.dset_edit.meta["parser"]["file_path"],
            )

        log.info(
            "Calculating satellite position/velocity (precise) based on SP3 precise orbit file {}.",
            ", ".join(self.dset_edit.meta["parser"]["file_path"]),
        )

        dset.vars["orbit"] = self.name
        dset.num_obs = len(self.time)
        dset.add_time("time", val=self.time, scale=self.time.scale)
        dset.add_text("satellite", val=self.satellite)
        dset.add_text("system", val=self.system)

        sat_pos = np.zeros((dset.num_obs, 3))
        sat_vel = np.zeros((dset.num_obs, 3))
        ref_time = dset.time[0]  # Reference epoch used for interpolation

        # Loop over all given satellites
        for sat in set(self.satellite):

            log.debug("Interpolation for satellite: {}", sat)

            # Get array with information about, when observation are available for the given satellite (indicated by
            # True)
            idx = dset.filter(satellite=sat)
            orb_idx = self.dset_edit.filter(satellite=sat)

            # Interpolation for given observation epochs (transmission time)
            sat_pos[idx], sat_vel[idx] = interpolation.interpolate_with_derivative(
                self.dset_edit.time[orb_idx].gps.sec_to_reference(ref_time),
                self.dset_edit.sat_pos.itrs[orb_idx],
                dset.time[idx].gps.sec_to_reference(ref_time),
                kind="lagrange",
                window=10,
                dx=0.5,
            )

            if np.isnan(np.sum(sat_pos[idx])) or np.isnan(np.sum(sat_vel[idx])):
                log.fatal(
                    "NaN occurred by determination of precise satellite position and velocity for satellite {}.", sat
                )

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
        #                     ' '.join([str(v) for v in dset.sat_posvel.itrs_pos[num_obs]]),
        #                     dset.gnss_satellite_clock[num_obs],
        #                     dset.gnss_relativistic_clock[num_obs])
        # -DEBUG

        return dset

    def satellite_clock_correction(self):
        """Determine satellite clock correction based on precise satellite clock product

        The GNSS satellite clock bias is read from RINEX clock files. Afterwards the satellite clock bias is determined
        via a cubic interpolation for the observation time.

        TODO:
            * Beware of the extrapolation (bounds_error=False in interpolate).
            * Check for satellite clock interpolation in:
              "Hesselbarth, A.: Statische und kinematische GNSS-Auswertung mittels PPP, 2011"

        Returns:
            numpy.ndarray:    GNSS satellite clock corrections for each observation
        """
        correction = np.zeros(self.dset.num_obs)
        sat_transmission_time = self.dset.time.gps.gpssec

        # Get precise GNSS satellite clock values
        clock_product = config.tech.get("clock_product", default="clk").str
        if clock_product == "sp3":
            all_sat_clk = self.dset_edit
        elif clock_product == "clk":
            all_sat_clk = data.Dataset(
                rundate=self.dset.rundate, tech=None, stage=None, dataset_name="gnss_sat_clk", dataset_id=0, empty=True
            )
            parser = parsers.parse("rinex_clk", rundate=self.dset.rundate)
            parser.write_to_dataset(all_sat_clk)  # TODO Read RINEX clock file, from day before and day after.
            #     Needed for interpolation. Add handling if these clk-files
            #     are not available. Remove first and last observations?
            #     If precise clock products are not available broadcast
            #     ephemeris should be used.
        else:
            log.fatal(
                "Unknown clock product '{}'. Configuration option 'clock_product' can only be 'sp3' or 'clk'.",
                clock_product,
            )

        # Loop over all satellites given in configuration file
        for sat in self.dset.unique("satellite"):

            # Skip satellites, which are not given in RINEX clock file
            if sat not in all_sat_clk.unique("satellite"):
                # TODO: Maybe satellite is available in SP3 file, which includes also
                #      satellite clock bias, but only for every 15 min instead of
                #      every 5 min (or 30 s by use of igs<wwwwd>.clk_30s).
                continue

            idx = self.dset.filter(satellite=sat)
            clk_idx = all_sat_clk.filter(satellite=sat)

            # Interpolation of GNSS precise satellite clock values
            # TODO: Check if interpolation method is ok.
            sat_clock_bias_ip = interpolate.interp1d(
                all_sat_clk.time.gps.gpssec[clk_idx],
                all_sat_clk.sat_clock_bias[clk_idx],
                axis=0,
                kind="cubic",
                bounds_error=False,
                fill_value=all_sat_clk.sat_clock_bias[clk_idx][-1],
            )
            correction[idx] = sat_clock_bias_ip(sat_transmission_time[idx])

        return correction
