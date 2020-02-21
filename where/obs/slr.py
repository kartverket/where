"""Populates a SLR Dataset with information from observation files and a priori sources

"""
from datetime import datetime, timedelta
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math import interpolation
from midgard.dev import exceptions
from midgard.math.unit import Unit

# Where imports
from where import apriori
from where import parsers
from where.ext import sofa
from where.lib import config
from where.lib import log

TECH = __name__.split(".")[-1]


@plugins.register
def write_to_dataset(dset, rundate=None, obs_format=None, **obs_args):

    obs_format = config.tech.get("obs_format", section=TECH, value=obs_format).str
    log.info(f"Reading observation file in {obs_format} format")

    file_vars1 = config.create_file_vars(rundate, TECH, **obs_args)
    last_date_to_read = rundate + timedelta(days=config.tech.arc_length.float + 1)
    parser1 = parsers.parse_key(f"slr_obs_{obs_format}", file_vars1)
    file_vars2 = config.create_file_vars(last_date_to_read, TECH, **obs_args)
    parser2 = parsers.parse_key(f"slr_obs_{obs_format}", file_vars2)

    if parser1.data_available and parser2.data_available:
        data = _write_to_dataset(parser1, parser2, dset, rundate)
        _write_met_to_dataset(dset, data, rundate)
    elif parser2.data_available and not parser2.data_available:
        raise exceptions.MissingDataError(
            f"No observation file in {obs_format} format found for {last_date_to_read.month}"
        )
    else:
        raise exceptions.MissingDataError(f"No observation file in {obs_format} format found for {rundate}")


def _write_to_dataset(parser1, parser2, dset, rundate):
    """Store SLR data in a dataset"""

    data_all1 = parser1.as_dict()
    data_all2 = parser2.as_dict()
    if parser1.file_path == parser2.file_path:
        collection = [data_all1]
    else:
        collection = [data_all1, data_all2]

    # Meta information
    dset.meta["tech"] = "slr"
    dset.meta.add("file", parser1.file_path.stem, section="input")
    dset.meta.add("file", parser2.file_path.stem, section="input")
    dset.meta.add("type", config.tech.obs_format.str.upper(), section="input")

    # Make new dict "obs_data" containing only data in relevant time interval:
    arc_length = config.tech.arc_length.float
    rundate_datetime = datetime(rundate.year, rundate.month, rundate.day)
    obs_data = dict()
    for data_all in collection:
        for i, x in enumerate(data_all["meta"]["obs_time"]):
            if rundate_datetime <= x < rundate_datetime + timedelta(days=arc_length):
                for key in ("meta", "obs", "obs_str"):
                    for field, val in data_all[key].items():
                        obs_data.setdefault(key, dict()).setdefault(field, list()).append(val[i])

        data_all.pop("meta")
        data_all.pop("obs")
        data_all.pop("obs_str")

        for key in data_all.keys():
            if key.startswith("met_"):
                for key2, val in data_all[key].items():
                    obs_data.setdefault(key, dict()).setdefault(key2, list()).append(val)
            elif key.startswith("satellite_"):
                # TODO: Use this information in the future?
                continue
            elif key.startswith("station_"):
                # TODO: Use this information in the future?
                continue
            else:
                log.fatal(f"Unknown data type{key}")

    obs_date = obs_data["meta"]["obs_date"]
    time = [obs_date[i] + timedelta(seconds=obs_data["meta"]["obs_sec"][i]) for i in range(0, len(obs_date))]
    dset.num_obs = len(obs_data["meta"]["obs_time"])
    dset.add_time("time", val=time, scale="utc", fmt="datetime")
    dset.add_text(val=obs_data["meta"]["station"], name="station")
    dset.add_text(val=obs_data["meta"]["satellite"], name="satellite")
    dset.add_float(val=obs_data["meta"]["bin_rms"], unit="picoseconds", name="bin_rms")
    # Positions
    trf = apriori.get("trf", time=dset.time)
    for station in dset.unique("station"):
        trf_site = trf[station]
        station_pos = trf_site.pos.trs.val
        log.debug(f"Station position for {station} ({trf_site.name}) is (x,y,z) = {station_pos.mean(axis=0)}")
        domes = trf_site.meta["domes"]
        obs_data["pos_" + station] = station_pos
        obs_data["station-other_" + station] = dict(domes=domes, cdp=station, site_id=station)
    dset.add_position(
        "site_pos",
        time=dset.time,
        system="trs",
        val=np.array([obs_data["pos_" + s][idx] for idx, s in enumerate(dset.station)]),
    )
    # Station data
    sta_fields = set().union(*[v.keys() for k, v in obs_data.items() if k.startswith("station_")])
    for field in sta_fields:
        dset.add_float(field, val=np.array([float(obs_data["station_" + s][field]) for s in dset.station]))
    sta_fields = set().union(*[v.keys() for k, v in obs_data.items() if k.startswith("station-other_")])
    for field in sta_fields:
        dset.add_text(field, val=[obs_data["station-other_" + s][field] for s in dset.station])

    # Station meta
    station_keys = sorted([k for k, v in obs_data.items() if k.startswith("station-other_")])
    pos_keys = sorted([k for k, v in obs_data.items() if k.startswith("pos_")])

    for sta_key, pos_key in zip(station_keys, pos_keys):
        sta_name = sta_key.replace("station-other_", "")
        cdp = obs_data[sta_key]["cdp"]
        dset.meta.add(sta_name, "site_id", cdp)
        longitude, latitude, height, _ = sofa.iau_gc2gd(2, obs_data[pos_key][0, :])  # TODO: Reference ellipsoid
        dset.meta.add("cdp", cdp, section=sta_name)
        dset.meta.add("site_id", cdp, section=sta_name)
        dset.meta.add("domes", obs_data[sta_key]["domes"], section=sta_name)
        dset.meta.add("marker", " ", section=sta_name)
        dset.meta.add("description", " ", section=sta_name)
        dset.meta.add("longitude", longitude, section=sta_name)
        dset.meta.add("latitude", latitude, section=sta_name)
        dset.meta.add("height", height, section=sta_name)

    # Satellite data
    sat_fields = set().union(*[v.keys() for k, v in obs_data.items() if k.startswith("satellite_")])
    for field in sat_fields:
        dset.add_float(field, val=np.array([float(obs_data["satellite_" + s][field]) for s in dset.satellite]))

    # Observations
    # In the dataset, obs_time is seconds since rundate:
    v = [
        (obs_data["meta"]["obs_date"][i] - rundate_datetime).total_seconds() + obs_data["meta"]["obs_sec"][i]
        for i in range(0, dset.num_obs)
    ]

    obs_data["obs"].pop("obs_time")
    dset.add_float("obs_time", val=v)
    for field, values in obs_data["obs"].items():
        dset.add_float(field, val=np.array(values))

    for field, values in obs_data["obs_str"].items():
        dset.add_text(field, val=values)

    return obs_data


def _write_met_to_dataset(dset, data, rundate):
    """Write the meteorological data from the parser to the dataset
    """
    data = _interpolate_meteorological_data(dset, data, rundate)

    met_fields = set().union(*[v.keys() for k, v in data.items() if k.startswith("met_")])
    for field in met_fields:
        dset.add_float(field, val=np.diag([data["met_" + s][field] for s in dset.station]))


def _interpolate_meteorological_data(dset, data, rundate):
    """Calculate temperature, humidity and pressure at observation epochs

    Meteorological data are calculated at observation epochs by interpolating in the data given on the observation
    file for each station.

    Missing meteorological data are currently not handled.
    """
    rundate = datetime(rundate.year, rundate.month, rundate.day)
    for field, station in [(f, f[4:]) for f in data.keys() if f.startswith("met_")]:
        log.debug(f"Meteorological data available for station {station}")

        met_time = data[field].pop("met_time")
        flat_list = [item for sublist in met_time for item in sublist]
        met_time_float = np.array([(flat_list[i] - rundate).total_seconds() for i in range(0, len(flat_list))])
        met_time_unique, met_index = np.unique(met_time_float, return_index=True)

        diff = len(met_time_float) - len(met_time_unique)
        if diff > 0:
            log.dev(f"Removed duplicate met data for station {station}")
            log.dev("Do this for the actual obs data also!")
        if len(met_time_unique) == 1:
            for met_type in data[field].keys():
                data[field][met_type] = np.repeat(data[field][met_type][0], dset.num_obs)
            continue

        # Extrapolation one month before/after
        # (this is overkill, most of these values will be removed later when taking the diagonal)
        min_time = min(met_time_unique) - 31 * 86400
        max_time = max(met_time_unique) + 31 * 86400
        met_time_unique = np.hstack((np.array(min_time), met_time_unique, np.array(max_time)))

        for met_type in data[field].keys():
            met_data_array = data[field][met_type]
            flat_list = [item for sublist in met_data_array for item in sublist]
            met_data_array = np.array([flat_list[i] for i in met_index])
            met_data_array = np.hstack((met_data_array[0], met_data_array, met_data_array[-1]))
            data[field][met_type] = interpolation.interpolate(
                met_time_unique, met_data_array, dset.obs_time, kind="cubic"
            )

    return data
