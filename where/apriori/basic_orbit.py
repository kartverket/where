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

    orbit = apriori.get('basic_orbit', rundate=rundate, days_after=1, days_before=1)

After the orbit dictionary is created, you can get the different functions::

    orbit["G10"]["pos"](dset.time)
    orbit["G10"]["vel"](dset.time)
    
The positions and velocities are computed in a terrestrial reference system (trs).

The time input argument to the functions must be covered by the time period in the apriori.get call.
Otherwise the functions will raise a MissingDataError exception.

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
def get_orbit(rundate, file_key=None, days_before=1, days_after=1):
    """
    Returns a dictionary with interpolation functions for position and velocity for each
    GNSS satellite. Reads sp3 files as defined in the file_key gnss_orbit_sp3. The files
    for the dates [rundate - days_before, rundate + days_after] are read and used to create
    the interpolation functions.
    
    Args:
        rundate (date):     Date of model run.
        days_before (int):  Number of days of sp3 files to read before the rundate
        days_after (int):   Number of days of sp3 files to read after the rundate
    
    Returns:
        dict:               Dictionary with interpolation functions. Structure:
                            {sat_name: {"pos": pos_func, {"vel": vel_func}}}
                            pos_func and vel_func requires a Time object as input
    """
    file_key = "gnss_orbit_sp3" if file_key is None else file_key
    
    date_to_read = rundate - timedelta(days=days_before)
    orb_data = {}
    
    # Read the files for all the days and collect it in orb_data
    while date_to_read <= rundate + timedelta(days=days_after):
        file_vars=config.date_vars(date_to_read)

        parser = parsers.parse_key("gnss_orbit_sp3", file_vars=file_vars)
        if parser.data is None:
            log.warn(f"Missing data from {parser.file_path}")
        else:
            log.info(f"Parsed precise orbit file {parser.file_path}")
        
        satellite = np.array(parser.data["satellite"])
        satellites = np.unique(satellite)
        time = np.array(parser.data["time"]) # GPS time in isot format
        sat_pos = np.array(parser.data["sat_pos"])
        
        for sat in satellites:
            sat_dict = orb_data.setdefault(sat, {})
            sat_dict.setdefault("time", [])
            sat_dict.setdefault("xyz", [])
            sat_dict["system"] = sat[0]
            idx = satellite == sat
            # Remove last epoch of the day because it is also included in the start of the next day
            sat_dict["time"] += time[idx][:-1].tolist()
            sat_dict["xyz"] += sat_pos[idx][:-1].tolist()
    
        date_to_read += timedelta(days=1)
    
    
    # Create interpolation functions for position and velocity for each satellite    
    orb = {}
    satellites = list(orb_data.keys())
    for sat in satellites:
        sat_dict = orb.setdefault(sat, {})
        sat_time = Time(orb_data[sat]["time"], fmt="isot", scale="gps")
        test = _time_to_jd2(sat_time)
        sat_dict["pos"] = _get_position_func(sat_time, orb_data[sat]["xyz"])
        sat_dict["vel"] = _get_velocity_func(sat_dict["pos"]) 
    
    return orb  

def _time_to_jd2(time):
    """Convert time object to fraction of day relative to first epoch.
    
    Orbits are given in GPS-time.
    """
    ref_time = time.gps[0].jd1
    ref_jd1 = time.gps.jd1 - ref_time
    jd2 = time.gps.jd2
    return ref_jd1 + jd2

def _get_position_func(sat_time, sat_pos):
    def position(time):
        """time: where.data.time.Time"""
        # Lagrange interpolation with 10 points: https://gssc.esa.int/navipedia/index.php/Precise_GNSS_Satellite_Coordinates_Computation
        func = lagrange(_time_to_jd2(sat_time), np.array(sat_pos), window=10)
        try:
            result = func(_time_to_jd2(time))
        except ValueError:
            raise exceptions.MissingDataError(f"Some orbit data missing for {sat_time[0]}-{sat_time[-1]}")
        return result
    return position

def _get_velocity_func(pos_func):
    def velocity(time):
        """time: where.data.time.Time"""
        # Estimate velocity based on position right before and after given epoch
        dt = TimeDelta(np.array([1/86400]*len(time)), scale="gps", fmt="days") # 1 second
        return (pos_func(time.gps - dt) - pos_func(time.gps + dt))/(2 * dt.val[:, None])
    return velocity
