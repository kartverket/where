"""Provides interpolated data from 

Description:
------------

from where import apriori
ntapl = apriori.get("non_tidal_atmospheric_loading_yearly_station", time=time)
dup = ntapl[station]["up"](time)
"""

import numpy as np
from scipy.interpolate import interp1d
from datetime import timedelta

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where import parsers


@plugins.register
def get_non_tidal_atmospheric_loading_yearly_station(time):
    """Read station data files relevant for the given time epochs

    Args:
        time (Time):    observation epochs

    Returns:
        A function that can interpolate
    """
    start = time.utc.datetime.date() if len(time) == 1 else min(time.utc.datetime).date()
    end = time.utc.datetime.date() if len(time) == 1 else max(time.utc.datetime).date()
    
    date_to_read = start
    loading_data = dict()

    while date_to_read <= end:
        file_vars = dict(config.date_vars(date_to_read))
        parser = parsers.parse_key(file_key="non_tidal_atmospheric_loading_yearly_station", file_vars=file_vars)
        data_chunk = parser.as_dict()
        for sta, sta_data in data_chunk.items():
            loading_data.setdefault(sta, {})
            for k in sta_data.keys():
                loading_data[sta][k] = loading_data[sta].get(k, []) + data_chunk[sta][k]
        date_to_read += timedelta(days=1)
    
    funcs = dict()
    for sta, sta_data in loading_data.items():
        funcs.setdefault(sta, {})
        mjd = sta_data.pop("mjd")
        # Use linear interpolation between each datapoint
        funcs[sta] = {k: interp1d(mjd,v, fill_value="extrapolate") for k, v in sta_data.items()}
    return funcs


@plugins.register
def get_non_tidal_atmospheric_loading(time):
    """Read gridded data files relevant for the given time epochs

    Args:
        time (Time):    observation epochs

    Returns:
        A function that can interpolate
    """
    min_time = min(time.utc).datetime
    max_time = max(time.utc).datetime
    start_hour = 6 * (min_time.hour // 6)
    start = min_time.replace(hour=start_hour, minute=0, second=0, microsecond=0)
    end_hour = 6 * (max_time.hour // 6)
    end = max_time.replace(hour=end_hour, minute=0, second=0, microsecond=0) + timedelta(hours=6)

    dt_to_read = start
    data = dict(up={}, east={}, north={})

    while dt_to_read <= end:
        file_vars = dict(config.date_vars(dt_to_read))

        data_chunk = parsers.parse_key(file_key="non_tidal_atmospheric_loading", file_vars=file_vars).as_dict()
        if data_chunk:
            data["up"][dt_to_read] = (data_chunk["lat"], data_chunk["lon"], data_chunk["up"])
            data["east"][dt_to_read] = (data_chunk["lat"], data_chunk["lon"], data_chunk["east"])
            data["north"][dt_to_read] = (data_chunk["lat"], data_chunk["lon"], data_chunk["north"])
        dt_to_read += timedelta(hours=6)

    return {k: create_interpolator(v) for k, v in data.items()}


def create_interpolator(data):
    """Creates interpolator for gridded dataset

    The interpolation between the in space is done with an bivariate spline interpolation.
    There is not an interpolation in time. Data from the last valid epoch is used.

    Args:
        data:    Dictionary with griddded values at given latitudes and longitudes.
                 The values are valid for 6 hours

                       data = { '<datetime>': (<lat list>, <lon list>, <value array (# of lat x # of lon)>),
                                     '<dateime>': (<lat list>, <lon list>, <value array (# of lat x # of lon)>),
                                     ...
                                    }

    @todo: Why is a bivariate spline interpolator used instead of bilinear interpolation?
           Explanation? Is there a significant difference to the bilinear interpolation?
           Maybe test needed scipy.interpolate.RectBivariateSpline against
           bilinear interpolation (Is that scipy.interpolate.interp2d?)?

    Returns:
        Interpolator function
    """
    # Create interpolator for gridded data
    interpolator = {rundate: RectBivariateSpline(lon, lat, values.T) for rundate, (lat, lon, values) in data.items()}

    def interpolate(time, longitude, latitude):
        """Interpolates in space and time for gridded data

        The interpolator function supports both multiple and single values (input can be either array/list or
        float/int).

        Args:
            dt_utc:     Datetime object(s) in UTC.
            longitude:  Longitude(s) in [rad]
            latitude:   Latitude(s) in [rad]

        Returns:
            Interpolated value(s)
        """
        if isinstance(longitude, (np.ndarray, list)):
            values = list()
            for obstime, lon, lat in zip(time.utc.datetime, longitude, latitude):
                starttime = obstime.replace(hour=6 * (obstime.hour // 6), minute=0, second=0, microsecond=0)
                # endtime = starttime + timedelta(hours=6)
                # fraction = (obstime - starttime).total_seconds() / 21600  # 21600 seconds per 6 hours

                # Linear time interpolation between hourly  files
                # value = interpolator[starttime](lon, lat, grid=False) + fraction * (
                #    interpolator[endtime](lon, lat, grid=False) - interpolator[starttime](lon, lat, grid=False)
                # )
                # value =
                values.append(interpolator[starttime](lon, lat, grid=False))
            return np.array(values)
        elif isinstance(longitude, (float, int)):
            dt_utc = time.utc.datetime
            starttime = dt_utc.replace(hour=6 * (dt_utc.hour // 6), minute=0, second=0, microsecond=0)
            # endtime = starttime + timedelta(hours=6)
            # fraction = (dt_utc - starttime).total_seconds() / 21600  # 21600 seconds per 6 hours

            # Linear time interpolation between hourly files
            # value = interpolator[starttime](longitude, latitude, grid=False) + fraction * (
            #    interpolator[endtime](longitude, latitude, grid=False)
            #    - interpolator[starttime](longitude, latitude, grid=False)
            # )
            # value =
            return interpolator[starttime](longitude, latitude, grid=False)
        else:
            log.fatal(f"Input {type(longitude)} is not a list, array, float or int")

    return interpolate
