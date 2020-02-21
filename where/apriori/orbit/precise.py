"""Handling of GNSS precise orbits

Description:
------------
The module includes a class for handling apriori GNSS precise orbits.

Example:
    from where import apriori

    # Get broadcast ephemeris object
    precise = apriori.get('orbit', apriori_orbit='precise')

    # Write calculated Dataset to file
    precise.dset.write()



"""
# Standard library imports
from datetime import timedelta
from typing import List, Tuple

# External library imports
from midgard.math import interpolation
import numpy as np
import pandas as pd
from scipy import interpolate

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import cleaners
from where.apriori import orbit
from where.data import dataset3 as dataset
from where.lib import config
from where.lib import log
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

    Methods:
        satellite_clock_correction():  Determine satellite clock correction
        _calculate():           Calculate precise orbit and satellite clock correction for given observation epochs
        _edit():                Edit precise orbit file data and save it in a Dataset
        _read():                Read precise orbit file data and save it in a Dataset

    """

    name = "precise"

    def __init__(self, rundate, file_key=None, file_path=None, day_offset=1, **kwargs):
        """Set up a new PreciseOrbit object, does not parse any data

        TODO: Remove dependency on rundate, use time to read correct files. (What to do with dataset?)

        Args:
            rundate (date):     Date of model run.
            file_key (str):     Key to the precise orbit file that will be read.
            file_path (pathlib.PosixPath):  File path to SP3 orbit file.
            day_offset (int):   Day offset used to calculate the number of days to read.
        """
        super().__init__(rundate=rundate)
        self.file_key = "gnss_orbit_sp3" if file_key is None else file_key
        self.file_path = file_path
        self.day_offset = day_offset

    def _read(self, dset_raw):
        """Read SP3 orbit file data and save it in a Dataset

        In addition to the given date, we read data for the day before and after. This is needed to carry out correct
        orbit interpolation at the start and end of a day.

        TODO:
        How well fits the orbits from day to day? Is it necessary to align the orbits?

        Args:
            dset_raw (Dataset):   Dataset representing raw data from apriori orbit files
        """
        date_to_read = dset_raw.analysis["rundate"] - timedelta(days=self.day_offset)
        file_paths = list()

        # Loop over days to read
        while date_to_read <= dset_raw.analysis["rundate"] + timedelta(days=self.day_offset):
            if self.file_path is None:
                file_path = config.files.path(self.file_key, file_vars=config.date_vars(date_to_read))
            else:
                file_path = self.file_path

            log.debug(f"Parse precise orbit file {file_path}")

            # Generate temporary Dataset with orbit file data
            dset_temp = dataset.Dataset(rundate=date_to_read, pipeline=dset_raw.vars["pipeline"], stage="temporary")
            parser = parsers.parse(parser_name="orbit_sp3", file_path=file_path, rundate=date_to_read)
            parser.write_to_dataset(dset_temp)
            file_paths.append(str(parser.file_path))

            # Extend Dataset dset_raw with temporary Dataset
            date = date_to_read.strftime("%Y-%m-%d")
            dset_raw.update_from(dset_temp) if dset_raw.num_obs == 0 else dset_raw.extend(dset_temp, meta_key=date)
            dset_raw.meta.add("file_path", file_paths, section="parser")

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
                f"Remove following duplicated precise orbit entries: \n{precise_filtered[['satellite', 'time']][idx]}"
            )
        precise_filtered = precise_filtered.drop_duplicates(subset=["satellite", "time"], keep="last")

        # TODO hjegei: Possible future ...
        # dset_edit.copy_from(self.dset_raw)
        # dset_edit.reorder(precise_filtered.index.values)

        # Convert edited fields to Dataset
        precise_np = precise_filtered.values
        fields = precise_filtered.columns

        dset_edit.vars["orbit"] = self.name
        dset_edit.num_obs = precise_filtered.shape[0]
        dset_edit.meta.update(self.dset_raw.meta)

        for idx, field in enumerate(fields):
            if field in ["time"]:
                # Note: It is ensured time scale consistency by using time scale from self.dset_raw. This is important,
                #       because 'self.dset_raw.as_dataframe()' uses time scale, when self.dset_raw time object was
                #       initialized.
                dset_edit.add_time(field, val=precise_np[:, idx], scale=self.dset_raw.time.scale, fmt="datetime")
            elif field in ["satellite", "system"]:
                dset_edit.add_text(field, val=precise_np[:, idx])
            elif field in ["sat_clock_bias"]:
                dset_edit.add_float(field, val=precise_np[:, idx].astype(float))

        # Add satellite position values to Dataset
        sat_pos = np.stack(
            (
                precise_np[:, fields.get_loc("sat_pos_trs_x")],
                precise_np[:, fields.get_loc("sat_pos_trs_y")],
                precise_np[:, fields.get_loc("sat_pos_trs_z")],
            ),
            axis=1,
        )
        dset_edit.add_position("sat_pos", val=sat_pos, time=dset_edit.time, system="trs")

        return dset_edit

    def _calculate(self, dset_out: "Dataset", dset_in: "Dataset", time: str = "time") -> None:
        """Calculate precise orbits and satellite clock correction for given observation epochs

        As a first step observations are removed from unavailable satellites and for exceeding the interpolation
        boundaries. The input Dataset contains observation epochs for which the broadcast ephemeris and satellite
        clock correction should be determined. The satellite position is determined for each observation epoch by
        interpolating within the given SP3 orbit time entries. The satellite velocities are calculated based on
        satellite positions 0.5 second before and after the observation epoch.

        Args:
            dset_out (Dataset): Dataset representing calculated precise orbits with following fields:

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

            dset_in:  Input Dataset containing model data for which broadcast ephemeris should be determined.
            time: Define time fields to be used. It can be for example 'time' or 'sat_time'. 'time' is related to 
                  observation time and 'sat_time' to satellite transmission time.
        """
        # Clean orbits by removing unavailable satellites, unhealthy satellites and checking interpolation boundaries
        cleaners.apply_remover("gnss_clean_orbit", dset_in, orbit_flag="precise")

        # TODO: Another solution has to be found for satellites not given in SP3 file, e.g. use of broadcast
        #       ephemeris.

        log.info(
            f"Calculating satellite position/velocity (precise) based on SP3 precise orbit file "
            f"{', '.join(self.dset_edit.meta['parser']['file_path'])}"
        )

        sat_pos = np.zeros((dset_in.num_obs, 3))
        sat_vel = np.zeros((dset_in.num_obs, 3))
        ref_time = dset_in[time][0]  # Reference epoch used for interpolation

        # Loop over all given satellites
        for sat in set(dset_in.satellite):

            log.debug(f"Interpolation for satellite: {sat}")

            # Get array with information about, when observation are available for the given satellite (indicated by
            # True)
            idx = dset_in.filter(satellite=sat)
            orb_idx = self.dset_edit.filter(satellite=sat)

            if np.min(dset_in[time][idx].gps.mjd) < np.min(self.dset_edit.time[orb_idx].mjd):
                log.fatal(
                    f"Interpolation range is exceeded by satellite {sat} ({np.max(dset_in[time][idx].gps.datetime)} "
                    f"[epoch] < {np.max(self.dset_edit.time[orb_idx].gps.datetime)} [precise orbit])"
                )

            if np.max(dset_in[time][idx].gps.mjd) > np.max(self.dset_edit.time[orb_idx].mjd):
                log.fatal(
                    f"Interpolation range is exceeded by satellite {sat} ({np.max(dset_in[time][idx].gps.datetime)} "
                    f"[epoch] > {np.max(self.dset_edit.time[orb_idx].gps.datetime)} [precise orbit])"
                )

            # Interpolation for given observation epochs (transmission time)
            diff_time_points = ref_time.gps - self.dset_edit.time.gps[orb_idx]
            diff_time_obs = ref_time.gps - dset_in[time].gps[idx]
            sat_pos[idx], sat_vel[idx] = interpolation.interpolate_with_derivative(
                # self.dset_edit.time[orb_idx].gps.sec_to_reference(ref_time),
                diff_time_points.seconds,
                self.dset_edit.sat_pos.trs[orb_idx],
                # dset_in[time][idx].gps.sec_to_reference(ref_time),
                diff_time_obs.seconds,
                kind="lagrange",
                window=10,
                dx=0.5,
            )

            if np.isnan(np.sum(sat_pos[idx])) or np.isnan(np.sum(sat_vel[idx])):
                log.fatal(
                    f"NaN occurred by determination of precise satellite position and velocity for satellite {sat}"
                )

        # Copy fields from model data Dataset
        dset_out.num_obs = dset_in.num_obs
        dset_out.add_text("satellite", val=dset_in.satellite)
        dset_out.add_text("system", val=dset_in.system)
        dset_out.add_time("time", val=dset_in[time])
        dset_out.vars["orbit"] = self.name

        # Add float fields
        dset_out.add_float(
            "gnss_relativistic_clock", val=self.relativistic_clock_correction(sat_pos, sat_vel), unit="meter"
        )
        dset_out.add_float(
            "gnss_satellite_clock", val=self.satellite_clock_correction(dset_in, time=time), unit="meter"
        )

        # Add satellite position and velocity to Dataset
        dset_out.add_posvel("sat_posvel", time=dset_out.time, system="trs", val=np.hstack((sat_pos, sat_vel)))

        # +DEBUG
        # for num_obs  in range(0, dset_out.num_obs):
        #    print('DEBUG: ', dset_out.satellite[num_obs],
        #                     dset_out[time].gps.datetime[num_obs],
        #                     dset_out[time].gps.mjd_frac[num_obs]*24*3600,
        #                     ' '.join([str(v) for v in dset_out.sat_posvel.trs.pos[num_obs]]),
        #                     dset_out.gnss_satellite_clock[num_obs],
        #                     dset_out.gnss_relativistic_clock[num_obs])
        # -DEBUG

    def _get_nearest_sample_point(self, satellite: Tuple[str], time: "Time") -> List[int]:
        """Get nearest sample point of precise orbits for given observation epochs

        Args:
            satellite:  Array with satellites related to given observation epochs
            time:       Observation epochs given as Time object

        Returns:
            List with nearest precise orbit sample point indices for given observation epochs.
        """
        precise_idx = list()
        log.debug(f"Get nearest interpolation sample points for given precise orbits.")

        # Check if precise orbits are available
        not_available_sat = sorted(set(satellite) - set(self.dset_edit.satellite))
        if not_available_sat:
            log.fatal(
                f"The following satellites are not given in apriori precise orbit file "
                f"{', '.join(self.dset_edit.meta['parser']['file_path'])}: {', '.join(not_available_sat)}"
            )

        if time.size > 1:
            # Determine nearest precise orbit sample point for a given satellite and observation epoch
            for sat, t in zip(satellite, time):

                idx = self.dset_edit.filter(satellite=sat)
                diff = t.gps.mjd - self.dset_edit.time.gps.mjd[idx]
                nearest_idx = np.array([abs(diff)]).argmin()

                precise_idx.append(idx.nonzero()[0][nearest_idx])
        else:
            for sat in satellite:

                idx = self.dset_edit.filter(satellite=sat)
                diff = time.gps.mjd - self.dset_edit.time.gps.mjd[idx]
                nearest_idx = np.array([abs(diff)]).argmin()

                precise_idx.append(idx.nonzero()[0][nearest_idx])

        return precise_idx

    def satellite_clock_correction(self, dset: "Dataset", time: str = "time") -> np.ndarray:
        """Determine satellite clock correction based on precise satellite clock product

        The GNSS satellite clock bias is read from RINEX clock files. Afterwards the satellite clock bias is determined
        via a cubic interpolation for the observation time.

        TODO:
            * Beware of the extrapolation (bounds_error=False in interpolate).
            * Check for satellite clock interpolation in:
              "Hesselbarth, A.: Statische und kinematische GNSS-Auswertung mittels PPP, 2011"

        Args:
            dset: A Dataset containing model data.
            time: Define time fields to be used. It can be for example 'time' or 'sat_time'. 'time' is related to 
                  observation time and 'sat_time' to satellite transmission time.

        Returns:
            GNSS satellite clock corrections for each observation
        """
        correction = np.zeros(dset.num_obs)
        sat_transmission_time = dset[time].gps.gps_ws.seconds

        # Get precise GNSS satellite clock values
        clock_product = config.tech.get("clock_product", default="clk").str
        if clock_product == "sp3":
            all_sat_clk = self.dset_edit
        elif clock_product == "clk":

            # TODO: File path information has to be improved, because 3 consecutive days are read.
            log.info(
                f"Calculating satellite clock correction (precise) based on RINEX clock file "
                f"{config.files.path(file_key='gnss_rinex_clk')}"
            )

            all_sat_clk = dataset.Dataset(rundate=dset.analysis["rundate"])
            parser = parsers.parse("rinex_clk", rundate=dset.analysis["rundate"])
            parser.write_to_dataset(all_sat_clk)  # TODO Read RINEX clock file, from day before and day after.
            #     Needed for interpolation. Add handling if these clk-files
            #     are not available. Remove first and last observations?
            #     If precise clock products are not available broadcast
            #     ephemeris should be used.
        else:
            log.fatal(
                f"Unknown clock product {clock_product!r}. "
                "Configuration option 'clock_product' can only be 'sp3' or 'clk'"
            )

        # Loop over all satellites given in configuration file
        for sat in dset.unique("satellite"):

            # Skip satellites, which are not given in RINEX clock file
            if sat not in all_sat_clk.unique("satellite"):
                # TODO: Maybe satellite is available in SP3 file, which includes also
                #      satellite clock bias, but only for every 15 min instead of
                #      every 5 min (or 30 s by use of igs<wwwwd>.clk_30s).
                continue

            idx = dset.filter(satellite=sat)
            clk_idx = all_sat_clk.filter(satellite=sat)

            # Interpolation of GNSS precise satellite clock values
            # TODO: Check if interpolation method is ok.
            sat_clock_bias_ip = interpolate.interp1d(
                all_sat_clk.time.gps.gps_ws.seconds[clk_idx],
                all_sat_clk.sat_clock_bias[clk_idx],
                axis=0,
                kind="cubic",
                bounds_error=False,
                fill_value=all_sat_clk.sat_clock_bias[clk_idx][-1],
            )
            correction[idx] = sat_clock_bias_ip(sat_transmission_time[idx])

        return correction
