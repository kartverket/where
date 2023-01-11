"""Populates a VLBI Dataset with information from observation files and a priori sources

"""
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.constant import constant
from midgard.math import ellipsoid

# where imports
from where import apriori
from where import parsers
from where.lib import log
from where.lib import config
from where.lib import exceptions
from where.ext import sofa

pipeline = __name__.split(".")[-1]


@plugins.register
def write_to_dataset(dset, rundate=None, obs_format=None, **obs_args):
    obs_format = config.tech.get("obs_format", section=pipeline, value=obs_format).str
    log.info(f"Reading observation file in {obs_format} format")
    session_code = dset.vars["session_code"]
    file_vars = config.create_file_vars(rundate, pipeline, session_code=session_code, **obs_args)
    parser = parsers.parse_key(f"vlbi_obs_{obs_format}", file_vars)

    if parser.data_available:
        _write_to_dataset(parser, dset, rundate, session_code)
    else:
        raise exceptions.MissingDataError(f"No observation file in {obs_format} format found for {rundate}")


def _write_to_dataset(parser, dset, rundate, session_code):

    data = parser.as_dict()
    units = data.get("meta", {}).get("units", {})

    # Session meta
    dset.meta.add("tech", "vlbi")
    dset.meta.add("file", parser.file_path.stem, section="input")
    dset.meta.add("type", config.tech.obs_format.str.upper(), section="input")

    master = apriori.get("vlbi_master_schedule", rundate=rundate)
    master_data = master.get(session_code, {})
    session_type = master_data.get("session_type", "")

    dset.meta.add("session_code", session_code, section="input")
    dset.meta.add("where_session_type", master.where_session_type(session_code), section="input")
    dset.meta.add("session_type", session_type, section="input")

    log.info(f"Session code: {session_code}")

    # Convert source names to official IERS names
    source_names = apriori.get("vlbi_source_names")
    iers_source_names = [source_names[src]["iers_name"] if src in source_names else src for src in data["source"]]
    # Replace the characeter "." with the letters "dot" in source names because "." has a special meaning in where
    data["source"] = iers_source_names

    # Replace spaces in station names with underscores to match official IVS name
    data["station_1"] = np.char.replace(data["station_1"], " ", "_")
    data["station_2"] = np.char.replace(data["station_2"], " ", "_")

    dset.num_obs = len(data["time"])
    dset.add_time("time", val=data.pop("time"), scale="utc", fmt="isot", write_level="operational")

    # Source directions
    crf = apriori.get("crf", time=dset.time)
    ra = np.array([crf[s].pos.right_ascension if s in crf else 0 for s in data["source"]])
    dec = np.array([crf[s].pos.declination if s in crf else 0 for s in data["source"]])
    dset.add_direction("src_dir", ra=ra, dec=dec, time=dset.time, write_level="operational")
    # Replace the characeter "." with the letters "dot" in source names because "." has a special meaning in where
    data["source"] = [s.replace(".", "dot") for s in iers_source_names]


    for field, values in data.items():
        values = np.array(values)
        if values.dtype.kind in {"U", "S"}:
            multiplier = -1 if field.endswith("_1") else 1
            dset.add_text(field, val=values, multiplier=multiplier, write_level="operational")
        elif values.dtype.kind in {"f", "i"}:
            multiplier = -1 if field.endswith("_1") else 1
            unit = units.get(field, None)
            dset.add_float(field, val=values, multiplier=multiplier, write_level="operational", unit=unit)
        elif values.dtype.kind in {"O"}:
            continue
        else:
            log.warn(f"Unknown datatype {values.dtype} for field {field}")


    # Station information
    log.info(f"Found stations: {', '.join(dset.unique('station'))}")
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
            trf_site = trf.closest(named_site.pos)
            cdp = trf_site.key
            ignore_stations = config.tech.ignore_station.stations.list
            logger = log.info if site in ignore_stations else log.warn
            logger(f"Undefined station name {site}. Assuming station is {trf_site.name} to get a cdp number.")

        data["pos_" + site] = trf_site.pos.trs.val
        _site_pos = np.mean(data[f"pos_{site}"], axis=0)
        log.debug(f"Using position {_site_pos} for {site} from {trf_site.source}")

        ivsname = station_codes[cdp]["name"]
        data["sta_" + site] = dict(site_id=cdp, cdp=cdp, ivsname=ivsname)

    # Positions
    itrs_pos_1 = np.array([data["pos_" + s][i, :] for i, s in enumerate(data["station_1"])])
    itrs_vel_1 = np.zeros((dset.num_obs, 3))
    dset.add_posvel(
        "site_pos_1",
        val=np.concatenate((itrs_pos_1, itrs_vel_1), axis=1),
        ellipsoid=ellipsoid.get(config.tech.reference_ellipsoid.str.upper()),
        system="trs",
        time=dset.time,
        # other=dset.src_dir,
        write_level="operational",
    )

    itrs_pos_2 = np.array([data["pos_" + s][i, :] for i, s in enumerate(data["station_2"])])
    itrs_vel_2 = np.zeros((dset.num_obs, 3))
    dset.add_posvel(
        "site_pos_2",
        val=np.concatenate((itrs_pos_2, itrs_vel_2), axis=1),
        ellipsoid=ellipsoid.get(config.tech.reference_ellipsoid.str.upper()),
        system="trs",
        time=dset.time,
        # other=dset.src_dir,
        write_level="operational",
    )

    # Compute aberrated source directions
    def aberrated_src_dir(site_pos):
        """See IERS2010 Conventions, equation 11.15"""
        site_vel_gcrs = site_pos.gcrs.vel.val
        eph = apriori.get("ephemerides", time=dset.time)
        vel = eph.vel_bcrs("earth") + site_vel_gcrs
        return (
            dset.src_dir.unit_vector
            + vel / constant.c
            - dset.src_dir.unit_vector * (dset.src_dir.unit_vector[:, None, :] @ vel[:, :, None])[:, :, 0] / constant.c
        )

    k_1 = aberrated_src_dir(dset.site_pos_1)
    dset.add_direction("abr_src_dir_1", val=k_1, system="gcrs", time=dset.time)
    dset.site_pos_1.other = dset.abr_src_dir_1

    k_2 = aberrated_src_dir(dset.site_pos_2)
    dset.add_direction("abr_src_dir_2", val=k_2, system="gcrs", time=dset.time)
    dset.site_pos_2.other = dset.abr_src_dir_2

    # Station data
    sta_fields = set().union(*[v.keys() for k, v in data.items() if k.startswith("sta_")])
    for field in sta_fields:
        dset.add_text(
            field + "_1", val=[data["sta_" + s][field] for s in data["station_1"]], multiplier=-1
        )  # write_level='analysis')
        dset.add_text(
            field + "_2", val=[data["sta_" + s][field] for s in data["station_2"]], multiplier=1
        )  # write_level='analysis')

    # Station meta
    station_keys = sorted([k for k, v in data.items() if k.startswith("sta_")])
    pos_keys = sorted([k for k, v in data.items() if k.startswith("pos_")])
    dset.meta["station"] = {}
    for sta_key, pos_key in zip(station_keys, pos_keys):
        sta_name = sta_key.replace("sta_", "")
        cdp = data[sta_key]["cdp"]
        ivsname = station_codes[cdp]["name"]
        longitude, latitude, height, _ = sofa.iau_gc2gd(2, data[pos_key][0, :])  # TODO: Reference ellipsoid

        dset.meta["station"].setdefault(ivsname, {})["cdp"] = cdp
        dset.meta["station"].setdefault(ivsname, {})["site_id"] = cdp
        dset.meta["station"].setdefault(ivsname, {})["domes"] = station_codes[cdp]["domes"]
        dset.meta["station"].setdefault(ivsname, {})["marker"] = station_codes[cdp]["marker"]
        dset.meta["station"].setdefault(ivsname, {})["description"] = station_codes[cdp]["description"]
        dset.meta["station"].setdefault(ivsname, {})["longitude"] = longitude
        dset.meta["station"].setdefault(ivsname, {})["latitude"] = latitude
        dset.meta["station"].setdefault(ivsname, {})["height"] = height
        if sta_name != ivsname:
            dset.meta["station"].setdefault(sta_name, {})["cdp"] = cdp
            dset.meta["station"].setdefault(sta_name, {})["site_id"] = cdp
            dset.meta["station"].setdefault(sta_name, {})["domes"] = station_codes[cdp]["domes"]
            dset.meta["station"].setdefault(sta_name, {})["marker"] = station_codes[cdp]["marker"]
            dset.meta["station"].setdefault(sta_name, {})["description"] = station_codes[cdp]["description"]
            dset.meta["station"].setdefault(sta_name, {})["longitude"] = longitude
            dset.meta["station"].setdefault(sta_name, {})["latitude"] = latitude
            dset.meta["station"].setdefault(sta_name, {})["height"] = height

    # Final cleanup
    # If there are more than 300 sources in a NGS-file the source names are gibberish
    bad_source_idx = ra == 0
    bad_sources = np.array(dset.source)[bad_source_idx]
    for s in np.unique(bad_sources):
        log.warn(f"Unknown source {s}. Observations with this source is discarded")
    dset.subset(np.logical_not(bad_source_idx))
