"""Populates a VLBI Dataset with information from observation files and a priori sources

"""
# Get info about session and date
# Find info about obs format
# Find obs_version
# Call correct parser
# Construct dataset
import re
import numpy as np

from where import apriori
from where import parsers
from where.lib import log
from where.lib import plugins
from where.lib import config
from where.lib import exceptions
from where.ext import sofa

TECH = __name__.split(".")[-1]


@plugins.register
def write_to_dataset(dset, rundate=None, session=None, obs_format=None, **obs_args):
    obs_format = config.tech.get("obs_format", section=TECH, value=obs_format).str
    log.info(f"Reading observation file in {obs_format} format")

    file_vars = config.create_file_vars(rundate, TECH, session=session, **obs_args)
    parser = parsers.parse_key(f"vlbi_obs_{obs_format}", file_vars)

    if parser.data_available:
        _write_to_dataset(parser, dset, rundate, session)
    else:
        raise exceptions.MissingDataError(f"No observation file in {obs_format} format found for {rundate}")


def _write_to_dataset(parser, dset, rundate, session):

    data = parser.as_dict()
    # TODO: units on fields

    # Convert source names to official IERS names
    source_names = apriori.get("vlbi_source_names")
    iers_source_names = [source_names[src]["iers_name"] if src in source_names else src for src in data["source"]]
    data["source"] = iers_source_names

    # Replace spaces in station names with underscores to match official IVS name
    data["station_1"] = np.char.replace(data["station_1"], " ", "_")
    data["station_2"] = np.char.replace(data["station_2"], " ", "_")

    dset.num_obs = len(data["time"])
    dset.add_time("time", val=data.pop("time"), scale="utc", format="isot", write_level="operational")
    for field, values in data.items():
        values = np.array(values)
        if values.dtype.kind in {"U", "S"}:
            dset.add_text(field, val=values, write_level="operational")
        elif values.dtype.kind in {"f", "i"}:
            dset.add_float(field, val=values, write_level="operational")
        elif values.dtype.kind in {"O"}:
            continue
        else:
            log.warn("Unknown datatype {} for field {}", values.dtype, field)

    # Source directions
    crf = apriori.get("crf", session=session)
    ra = np.array([crf[s].pos.crs[0] if s in crf else 0 for s in data["source"]])
    dec = np.array([crf[s].pos.crs[1] if s in crf else 0 for s in data["source"]])

    dset.add_direction("src_dir", ra=ra, dec=dec, write_level="operational")

    # Station information
    log.info("Found stations: {}", ", ".join(dset.unique("station")))
    trf = apriori.get("trf", time=dset.time)
    station_codes = apriori.get("vlbi_station_codes")
    dset.add_text(
        "baseline",
        val=np.array([f"{s1}/{s2}" for s1, s2 in zip(data["station_1"], data["station_2"])]),
        write_level="operational",
    )
    for site in dset.unique("station"):
        if site in station_codes:
            cdp = station_codes[site]["cdp"]
            trf_site = trf[cdp]
        else:
            named_site = trf.named_site(site)
            trf_site = trf.closest(named_site.pos, max_distance=5)
            cdp = trf_site.key
            ignore_stations = config.tech.ignore_station.stations.list
            if site in ignore_stations:
                log.info("Undefined station name {}. Assuming station is {}.", site, trf_site.name)
            else:
                log.warn("Undefined station name {}. Assuming station is {}.", site, trf_site.name)

        data["pos_" + site] = trf_site.pos.itrs
        log.debug("Using position {} for {} from {}", np.mean(data["pos_" + site], axis=0), site, trf_site.source)

        ivsname = station_codes[cdp]["name"]
        data["sta_" + site] = dict(site_id=cdp, cdp=cdp, ivsname=ivsname)

    # Positions
    itrs_pos_1 = np.array([data["pos_" + s][i, :] for i, s in enumerate(data["station_1"])])
    itrs_vel_1 = np.zeros((dset.num_obs, 3))
    dset.add_posvel(
        "site_pos_1",
        time="time",
        other="src_dir",
        itrs=np.concatenate((itrs_pos_1, itrs_vel_1), axis=1),
        write_level="operational",
    )
    itrs_pos_2 = np.array([data["pos_" + s][i, :] for i, s in enumerate(data["station_2"])])
    itrs_vel_2 = np.zeros((dset.num_obs, 3))
    dset.add_posvel(
        "site_pos_2",
        time="time",
        other="src_dir",
        itrs=np.concatenate((itrs_pos_2, itrs_vel_2), axis=1),
        write_level="operational",
    )

    # Station data
    sta_fields = set().union(*[v.keys() for k, v in data.items() if k.startswith("sta_")])
    for field in sta_fields:
        dset.add_text(
            field + "_1", val=[data["sta_" + s][field] for s in data["station_1"]]
        )  # write_level='analysis')
        dset.add_text(
            field + "_2", val=[data["sta_" + s][field] for s in data["station_2"]]
        )  # write_level='analysis')

    # Station meta
    station_keys = sorted([k for k, v in data.items() if k.startswith("sta_")])
    pos_keys = sorted([k for k, v in data.items() if k.startswith("pos_")])

    for sta_key, pos_key in zip(station_keys, pos_keys):
        sta_name = sta_key.replace("sta_", "")
        cdp = data[sta_key]["cdp"]
        ivsname = station_codes[cdp]["name"]
        longitude, latitude, height, _ = sofa.iau_gc2gd(2, data[pos_key][0, :])  # TODO: Reference ellipsoid
        dset.add_to_meta(ivsname, "cdp", cdp)
        dset.add_to_meta(ivsname, "site_id", cdp)
        dset.add_to_meta(ivsname, "domes", station_codes[cdp]["domes"])
        dset.add_to_meta(ivsname, "marker", station_codes[cdp]["marker"])
        dset.add_to_meta(ivsname, "description", station_codes[cdp]["description"])
        dset.add_to_meta(ivsname, "longitude", longitude)
        dset.add_to_meta(ivsname, "latitude", latitude)
        dset.add_to_meta(ivsname, "height", height)
        if sta_name != ivsname:
            dset.add_to_meta(sta_name, "cdp", cdp)
            dset.add_to_meta(sta_name, "site_id", cdp)
            dset.add_to_meta(sta_name, "domes", station_codes[cdp]["domes"])
            dset.add_to_meta(sta_name, "marker", station_codes[cdp]["marker"])
            dset.add_to_meta(sta_name, "description", station_codes[cdp]["description"])
            dset.add_to_meta(sta_name, "longitude", longitude)
            dset.add_to_meta(sta_name, "latitude", latitude)
            dset.add_to_meta(sta_name, "height", height)

    dset.meta["tech"] = "vlbi"
    dset.add_to_meta("input", "file", parser.file_path.stem)
    dset.add_to_meta("input", "type", config.tech.obs_format.str.upper())

    if "meta" not in data:
        master = apriori.get("vlbi_master_schedule", rundate=rundate)
        master_data = master.get((rundate.timetuple().tm_yday, session), {})
        dset.add_to_meta("input", "session_code", master_data.get("session_code", ""))
    else:
        dset.add_to_meta("input", "session_code", data["meta"].get("session_code", ""))

    reg_hits = re.search("\d", dset.meta["input"]["session_code"])
    num_idx = reg_hits.start() if reg_hits else len(dset.meta["input"]["session_code"])
    dset.add_to_meta("input", "session_type", dset.meta["input"]["session_code"][:num_idx])

    # Final cleanup
    # If there are more than 300 sources in a NGS-file the source names are gibberish
    bad_source_idx = ra == 0
    bad_sources = np.array(dset.source)[bad_source_idx]
    for s in np.unique(bad_sources):
        log.warn("Unknown source {}. Observations with this source is discarded", s)
    dset.subset(np.logical_not(bad_source_idx))
