"""Provides interpolated data from Vienna Mapping Function 1

Description:

Reads Vienna Mapping Function 1 coefficients, zenith hydrostatic and wet
delays and grid heights for data points and creates interpolation functions for each datatype.

The VMF1 gridded files are based on a global grid (2.0 x 2.5 degrees). The mapping
function coefficients "ah" and "aw" and the hydrostatic and wet zenith delay "zh"
and "zw" are given on a global grid with 2.0 degrees sampling from north to south
and 2.5 degrees sampling from west to east. For each parameter there are four
files per day, i.e. at 0, 6, 12, and 18 UT.

The file corresponding to the first 6 hours of the given rundate is read, and then
the file corresponding to the next 6 hours is read until all datafiles for the
given rundate and day before and after are read for each datatype ("ah", "aw", "zh",
"zw").

In addition the "orography_ell" file is read, which is including ellipsoidal
heights. All gridded VMF1 data refers to the ellipsoidal heights given in the file.

References:
-----------
    http://ggosatm.hg.tuwien.ac.at/DELAY/readme.txt

@todo which interpolators should be used?



"""

import numpy as np
from scipy.interpolate import RectBivariateSpline
from datetime import timedelta

# Where imports
from where.lib import config
from where.lib import plugins
from where import parsers

DATATYPE = ("ah", "aw", "zh", "zw")


@plugins.register
def get_vmf1(time):
    """Read VMF1 gridded data files relevant for the given time epochs

    Args:
        time (Time):    observation epochs

    Returns:
        A dictionary of functions that can interpolate in the VMF1 dataset.
    """
    data = dict()
    min_time = min(time.utc).datetime
    max_time = max(time.utc).datetime
    start_hour = 6 * (min_time.hour // 6)
    start = min_time.replace(hour=start_hour, minute=0, second=0, microsecond=0)
    end_hour = 6 * (max_time.hour // 6)
    end = max_time.replace(hour=end_hour, minute=0, second=0, microsecond=0) + timedelta(hours=6)

    for datatype in DATATYPE:
        dt_to_read = start
        vmf1_data = dict()

        while dt_to_read <= end:
            file_vars = dict(config.date_vars(dt_to_read), type=datatype)
            data_chunk = parsers.parse_key(file_key="vmf1", file_vars=file_vars).as_dict()
            vmf1_data[dt_to_read] = (data_chunk["lat"], data_chunk["lon"], data_chunk["values"])
            dt_to_read += timedelta(hours=6)
        data[datatype] = vmf1_data

    funcs = {k: vmf1_interpolator(v) for k, v in data.items()}

    data = parsers.parse_key(file_key="orography_ell").as_dict()
    funcs["ell"] = RectBivariateSpline(data["lon"], data["lat"], data["values"].T)
    return funcs


def vmf1_interpolator(vmf1_data):
    """Creates interpolator for VMF1 gridded dataset

    The interpolation between the VMF1 gridded datasets is done with an bivariate spline interpolation, whereas a
    linear interpolation in time is used between the 6 hourly VMF1 data files.

    Args:
        vmf1_data:    Dictionary with VMF1 gridded data given for the hours 0, 6, 12 and 18 UT (as datetime object)
                      over three consecutive days. Each hourly dictionary includes arrays with latitude, longitude and
                      the VMF1 gridded data.

                       vmf1_data = { '<1st hour>': (<lat list>, <lon list>, <value array (# of lat x # of lon)>),
                                     '<2nd hour>': (<lat list>, <lon list>, <value array (# of lat x # of lon)>),
                                     ...
                                   }

    @todo: Why is a bivariate spline interpolator used instead of bilinear interpolation?
           Explanation? Is there a significant difference to the bilinear interpolation?
           Maybe test needed scipy.interpolate.RectBivariateSpline against
           bilinear interpolation (Is that scipy.interpolate.interp2d?)?

    @todo: A linear time interpolation is used, this should be tested and maybe improved.
           For example VieVS uses Lagrange interpolation.

    Returns:
        Interpolator function
    """
    # Create interpolator for gridded VMF1 data
    interpolator = {
        rundate: RectBivariateSpline(lon, lat, values.T) for rundate, (lat, lon, values) in vmf1_data.items()
    }

    def interpolate(time, longitude, latitude):
        """Interpolates in space and time in the gridded vmf1 data

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
                endtime = starttime + timedelta(hours=6)
                fraction = (obstime - starttime).total_seconds() / 21600  # 21600 seconds per 6 hours

                # Linear time interpolation between hourly VMF1 files
                value = interpolator[starttime](lon, lat, grid=False) + fraction * (
                    interpolator[endtime](lon, lat, grid=False) - interpolator[starttime](lon, lat, grid=False)
                )
                values.append(value)
            return np.array(values)
        elif isinstance(longitude, (float, int)):
            dt_utc = time.utc.datetime
            starttime = dt_utc.replace(hour=6 * (dt_utc.hour // 6), minute=0, second=0, microsecond=0)
            endtime = starttime + timedelta(hours=6)
            fraction = (dt_utc - starttime).total_seconds() / 21600  # 21600 seconds per 6 hours

            # Linear time interpolation between hourly VMF1 files
            value = interpolator[starttime](longitude, latitude, grid=False) + fraction * (
                interpolator[endtime](longitude, latitude, grid=False)
                - interpolator[starttime](longitude, latitude, grid=False)
            )

            return value
        else:
            log.fatal(f"Input {type(longitude)} is not a list, array, float or int")

    return interpolate
