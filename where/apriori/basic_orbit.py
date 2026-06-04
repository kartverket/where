"""Get apriori GNSS orbits from sp3 files

Description:
------------

Reads GNSS orbits from the configured data source and creates interpolation functions for each satellite.
Only the position is used from the sp3 files. All other data in the orbit files are discarded.

The following functions are provided:

========== ============================================================
Function   Description
========== ============================================================
pos(Time)     Lagrange interpolation of position based on sp3 coordinates
vel(Time)     Estimation of velocity based on position
========== ============================================================

Example:
--------

To use the orbits, simply get them from the apriori-package specifying the days you need::

    orbit = apriori.get('basic_orbit', rundate=rundate, bounds_error=True, days_after=1, days_before=1)

After the orbit dictionary is created, you can get the different functions::

    orbit["G10"]["pos"](dset.time)
    orbit["G10"]["vel"](dset.time)

The input time argument must be a array-like.

The positions and velocities are computed in a terrestrial reference system (trs).

Unless bounds_error is set to False, the time input argument to the functions must be covered by the time period
in the apriori.get call. Otherwise the functions will raise a MissingDataError exception.

If a satellite is completely missing in the sp3 orbits for the given time period the call will result in a
normal KeyError.

"""

import numpy as np
from datetime import timedelta

# Midgard imports
from midgard.dev import plugins
from midgard.math.interpolation import lagrange

# Where imports
from where.data.time import Time, TimeDelta
from where.lib import config
from where.lib import exceptions
from where.lib import log
from where import parsers



@plugins.register
def get_orbit(rundate, file_key=None, bounds_error=False, days_before=1, days_after=1):
    """
    Returns a dictionary with interpolation functions for position and velocity for each
    GNSS satellite. Reads sp3 files as defined in the file_key gnss_orbit_sp3. The files
    for the dates [rundate - days_before, rundate + days_after] are read and used to create
    the interpolation functions.
    
    bounds_error is set to False by default. This is to avoid error message for boundary
    conditions when the orbit time scale is different from used time scale. This might
    happen because time.utc is not guaranteed to be identical to time.gps.utc due to the
    utc_tai conversion that is numerically symmertric.

    Args:
        rundate (date):     Date of model run.
        file_key:           Which file_key to read
        bounds_error: (bool)Flag to enable error message instead of extrapolation
        days_before (int):  Number of days of sp3 files to read before the rundate
        days_after (int):   Number of days of sp3 files to read after the rundate
    
    Returns:
        dict:               Dictionary with interpolation functions. Structure:
                            {sat_name: {"pos": pos_func, {"vel": vel_func}}}
                            pos_func and vel_func requires a Time object as input
    """
    file_key = "gnss_orbit_sp3" if file_key is None else file_key
    # TODO: Multiple file keys?
    return _orbit_from_sp3(rundate, file_key, bounds_error, days_before, days_after)


def _orbit_from_sp3(rundate, file_key, bounds_error, days_before, days_after):
    date_to_read = rundate - timedelta(days=days_before)
    orb_data = {}
    
    sat_vel = None
    parsed_files = []
    # Read the files for all the days and collect it in orb_data
    while date_to_read <= rundate + timedelta(days=days_after):
        file_vars=config.date_vars(date_to_read)

        parser = parsers.parse_key(file_key, file_vars=file_vars)
        if parser.data is None:
            log.warn(f"Missing data from {parser.file_path}")
            continue

        if parser.file_path not in parsed_files:
            log.info(f"Parsed precise orbit file {parser.file_path}")
            satellite = np.array(parser.data["satellite"])
            satellites = np.unique(satellite)
            time = np.array(parser.data["time"]) # time in isot format
            sat_pos = np.array(parser.data["sat_pos"])
            sat_vel = np.array(parser.data["sat_vel"]) if "sat_vel" in parser.data else None

            for sat in satellites:
                sat_dict = orb_data.setdefault(sat, {})
                sat_dict.setdefault("time", [])
                sat_dict.setdefault("xyz", [])
                sat_dict.setdefault("vxyz", [])
                sat_dict["system"] = sat[0]
                idx = satellite == sat
                # Remove last epoch of the day because it is also included in the start of the next day
                sat_dict["time"] += time[idx][:-1].tolist()
                sat_dict["xyz"] += sat_pos[idx][:-1].tolist()
                if sat_vel is not None:
                    sat_dict["vxyz"] += sat_vel[idx][:-1].tolist()

            # Avoid reading the same file twice if a single sp3 file contains data for more than one day
            parsed_files.append(parser.file_path)
        date_to_read += timedelta(days=1)


    # Create interpolation functions for position and velocity for each satellite
    orb = {}
    satellites = list(orb_data.keys())
    # Assume the time_sys in the last read file is representative for all read files
    time_sys = parser.meta["time_sys"]
    if time_sys not in ("GPS", "UTC"):
        log.warn(f"System Time Indicator {time_sys} in SP3 file {parser.file_path} is not supported.")
    time_sys = "gps" # Force time sys to be gps even though sp3 file says utc # lageos-1 debugging
    for sat in satellites:
        sat_dict = orb.setdefault(sat, {})
        sat_time = Time(orb_data[sat]["time"], fmt="isot", scale=time_sys.lower())
        sat_dict["pos"] = _get_position_func(sat_time, orb_data[sat]["xyz"], bounds_error=bounds_error)
        if sat_vel is not None:
            sat_dict["vel"] = _get_interpolated_velocity_func(sat_time, orb_data[sat]["vxyz"], bounds_error=bounds_error)
        else:
            sat_dict["vel"] = _get_derived_velocity_func(sat_dict["pos"])

    return orb

def _time_to_jd2(ref_time, time):
    """Convert time object to fraction of day relative to first epoch.

    Makes sure the time epochs use the same time scale as the ref_time

    Args:   ref_time: first time epoch of the satellite data used in interpolation
            time: epochs to convert to fraction of day
    """
    scale = ref_time.scale
    ref_jd1 = time.to_scale(scale).jd1 - ref_time.jd1
    jd2 = time.to_scale(scale).jd2
    output = ref_jd1 + jd2
    return ref_jd1 + jd2

def _get_position_func(sat_time, sat_pos, bounds_error):
    def position(time):
        """time: where.data.time.Time"""
        # Lagrange interpolation with 10 points: https://gssc.esa.int/navipedia/index.php/Precise_GNSS_Satellite_Coordinates_Computation
        ref_time = sat_time[0]
        func = lagrange(_time_to_jd2(ref_time, sat_time), np.array(sat_pos), bounds_error=bounds_error, window=10)
        try:
            result = func(_time_to_jd2(ref_time, time))
        except ValueError as err:
            raise exceptions.MissingDataError(f"Some orbit data missing for {sat_time[0]}-{sat_time[-1]} ({sat_time.scale}). {err}")
        return result
    return position

def _get_interpolated_velocity_func(sat_time, sat_vel, bounds_error):
    def velocity(time):
        """time: where.data.time.Time"""
        # Lagrange interpolation with 10 points: https://gssc.esa.int/navipedia/index.php/Precise_GNSS_Satellite_Coordinates_Computation
        ref_time = sat_time[0]
        func = lagrange(_time_to_jd2(ref_time, sat_time), np.array(sat_vel), bounds_error=bounds_error, window=10)
        try:
            result = func(_time_to_jd2(ref_time, time))
        except ValueError as err:
            raise exceptions.MissingDataError(f"Some orbit data missing for {sat_time[0]}-{sat_time[-1]}. {err}")
        return result
    return velocity

def _get_derived_velocity_func(pos_func):
    def velocity(time):
        """time: where.data.time.Time"""
        # Estimate velocity based on position right before and after given epoch
        dt = 1.0 # 1 second
        dt = TimeDelta(np.array([dt]*len(time)), scale=time.scale, fmt="seconds")
        return (pos_func(time - dt) - pos_func(time + dt))/(2 * dt.val[:, None])
    return velocity

