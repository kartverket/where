"""Provides interpolated data from Vienna Mapping Function 1

Description:

Usage:
from where import apriori

vmf1_station = apriori.get("vmf1_station", time=time)
ah = vmf1_station[station]["aw"](time.mjd)

References:
-----------
    http://ggosatm.hg.tuwien.ac.at/DELAY/readme.txt

"""
# External library imports
from scipy.interpolate import interp1d
from datetime import timedelta

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where import parsers


@plugins.register
def get_vmf1_station(time):
    """Read VMF1 station data files relevant for the given time epochs

    Args:
        time (Time):    observation epochs

    Returns:
        A dictionary of functions that can interpolate in the VMF1 dataset.
    """
    start = time.utc.datetime.date() if len(time) == 1 else min(time.utc.datetime).date()
    end = time.utc.datetime.date() if len(time) == 1 else max(time.utc.datetime).date()

    date_to_read = start
    vmf1_data = dict()

    while date_to_read <= end:
        file_vars = dict(config.date_vars(date_to_read))
        parser = parsers.parse_key(file_key="vmf1_station", file_vars=file_vars)
        data_chunk = parser.as_dict()
        for sta, sta_data in data_chunk.items():
            vmf1_data.setdefault(sta, {})
            for k in sta_data.keys():
                vmf1_data[sta][k] = vmf1_data[sta].get(k, []) + data_chunk[sta][k]
        date_to_read += timedelta(days=1)

    funcs = dict()
    for sta, sta_data in vmf1_data.items():
        funcs.setdefault(sta, {})
        mjd = sta_data.pop("mjd")
        # Use linear interpolation between each datapoint
        funcs[sta] = {k: interp1d(mjd,v, fill_value="extrapolate") for k, v in sta_data.items()}
    return funcs
