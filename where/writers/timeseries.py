"""Write key indicators from model run to timeseries dataset

Description:
------------

Todo

"""

# Standard library imports
from datetime import date, datetime
import itertools
import re
import sys

# External library imports
import numpy as np

# Where imports
from where.lib import config
from where import data
from where.lib import log
from where.lib import plugins


WRITER = __name__.split(".")[-1]


@plugins.register
def add_to_full_timeseries(dset):
    """Write some key variables to the full timeseries

    Args:
        dset:  Dataset, data for a model run.
    """
    dset_id = int(config.tech.timeseries.dataset_id.str.format(**dset.vars))
    dset_name = config.tech.timeseries.dataset_name.str.format(**dset.vars)
    dset_ts = data.Dataset(date(1970, 1, 1), dset.vars["tech"], "timeseries", dset_name, dset_id, session="")
    dset_session = data.Dataset.anonymous()

    # Add data to dset_session
    idx_fields = config.tech[WRITER].index.list
    field_values = [["all"] + dset.unique(f) for f in idx_fields]
    idx_values = dict(zip(idx_fields, zip(*itertools.product(*field_values))))
    # TODO: Remove combinations where filter leaves 0 observations

    num_obs = len(idx_values[idx_fields[0]])  # Length of any (in this case the first) field
    mean_epoch = dset.time.mean.utc
    rundate_str = dset.rundate.strftime(config.FMT_date)
    session = dset.dataset_name
    status = dset.meta.get("analysis_status", "unchecked")
    session_type = dset.meta.get("session_type", "")

    dset_session.num_obs = num_obs
    dset_session.add_time("time", val=[mean_epoch] * num_obs, scale="utc")
    dset_session.add_time("timestamp", val=[datetime.now()] * num_obs, scale="utc")
    dset_session.add_text("rundate", val=[rundate_str] * num_obs)
    dset_session.add_text("session", val=[session] * num_obs)
    dset_session.add_text("status", val=[status] * num_obs)
    dset_session.add_text("session_type", val=[session_type] * num_obs)

    for field, value in idx_values.items():
        dset_session.add_text(field, val=value)

    default_dset_str = f"{dset.vars['stage']}/{int(dset.dataset_id)}"
    dsets = {default_dset_str: dset}
    for method, cfg_entry in config.tech[WRITER].items():
        try:
            method_func = getattr(sys.modules[__name__], "method_{}".format(method))
        except AttributeError:
            log.warn("Method '{}' is unknown", method)
            continue

        for field_cfg in cfg_entry.as_list(split_re=", *"):
            field_out = re.sub("[ -/:]", "_", field_cfg)
            func, _, field_dset = field_cfg.rpartition(":")
            field_in, _, dset_str = field_dset.partition("-")
            func = func if func else field_in
            dset_str = dset_str if dset_str else default_dset_str
            if dset_str not in dsets:
                stage, _, dset_id = dset_str.partition("/")
                dset_id = int(dset_id) if dset_id else "last"
                dsets[dset_str] = data.Dataset(
                    dset.rundate,
                    tech=dset.vars["tech"],
                    stage=stage,
                    dataset_name=dset.dataset_name,
                    dataset_id=dset_id,
                )

            val, adder, unit = method_func(dsets[dset_str], field_in, idx_values, func)
            if adder:
                add_func = getattr(dset_session, adder)
                if val.ndim > 1:
                    add_func(field_out, val=val, shape=val.shape[1:], unit=unit)
                else:
                    add_func(field_out, val=val, unit=unit)

    # hack to get solved neq data into the time series:
    # TODO: unhack this :P Add as a method_neq instead?
    if "normal equation" in dset.meta:
        _add_solved_neq_fields(dset, dset_session, idx_values)

    # Filter timeseries dataset to remove any previous data for this rundate and session
    keep_idx = np.logical_not(dset_ts.filter(rundate=rundate_str, session=session))
    dset_ts.subset(keep_idx)

    # Extend timeseries dataset with dset_session and write to disk
    if dset_ts.num_obs:
        dset_ts.extend(dset_session)
    else:
        dset_ts.copy_from(dset_session)

    log.info("Updating timeseries dataset '{}'", dset_ts.description)
    dset_ts.write()


def _add_solved_neq_fields(dset, dset_session, idx_values):
    names = dset.meta["normal equation"]["names"]
    x = np.array(dset.meta["normal equation"]["solution"])
    units = np.array(dset.meta["normal equation"]["unit"])

    fields = np.unique([n.split("-")[0] for n in names])
    for field in fields:
        if "vlbi_src_dir" in field:
            # TODO
            continue
        idx_names = idx_values[list(idx_values.keys()).pop()]
        num_obs = len(idx_names)
        params = [f for f in names if f.startswith(field + "-")]
        param_units = [u for f, u in zip(names, units) if f.startswith(field + "-")]
        name2 = [p.split("-")[-1].split("_")[-1] for p in params]
        dim = len(np.unique(name2)) if not any([n2 in n for n2 in name2 for n in idx_names]) else 1
        val = np.full((num_obs, dim), np.nan, dtype=float)
        idx = np.array([any([n in p for p in params]) for n in idx_names], dtype=bool)
        mean = np.array([np.mean(state) for state, param in zip(x, names) if param in params]).reshape(-1, dim)
        if idx.any():
            val[idx] = mean
        else:
            val[0] = mean
        # all fields of the same type share the same unit, only the first entry of param_units is needed
        dset_session.add_float("neq_" + field, val=val, shape=val.shape[1:], unit=param_units[0])


#
# Methods, corresponds to keys in config.tech.timeseries
#
def method_index(dset, field, idx_values, func):
    return None, None, None


def method_func(dset, field, idx_values, func):
    calculate_func = getattr(sys.modules[__name__], func)
    return calculate_func(dset, field, idx_values)


def method_statistics(dset, field, idx_values, func):
    num_obs = len(idx_values[list(idx_values.keys()).pop()])
    val = np.full(num_obs, np.nan, dtype=float)
    val[0] = dset.meta["statistics"][field]

    return val, "add_float", None  # Todo: Can we add a unit here somehow?


def method_meta(dset, field, idx_values, func):
    # TODO: Geir Arne!
    return None, None, None


def method_text(dset, field, idx_values, func):
    text = np.array(["/".join(dset.unique(field, idx=i)) for i in _filter_each(dset, idx_values)])
    return text, "add_text", dset.unit(field)


def method_state(dset, field, idx_values, func):
    if "vlbi_src_dir" in field:
        # TODO
        return None, None, None
    idx_names = idx_values[list(idx_values.keys()).pop()]
    num_obs = len(idx_names)
    params = [f for f in dset.fields if f.startswith("state_" + field + "-") and not f.endswith("_sigma")]
    name2 = [p.split("-")[-1].split("_")[-1] for p in params]
    dim = len(np.unique(name2)) if not any([n2 in n for n2 in name2 for n in idx_names]) else 1
    val = np.full((num_obs, dim), np.nan, dtype=float)
    idx = np.array([any([n in p for p in params]) for n in idx_names], dtype=bool)
    mean = [dset.mean(p) for p in params]
    if mean:
        mean = np.array(mean).reshape(-1, dim)
    else:
        return None, None, None
    if idx.any():
        val[idx] = mean
    else:
        val[0] = mean

    # Get unit from config
    display_unit = config.tech[field].display_unit.str
    unit = config.tech[field].unit.str if not display_unit else display_unit

    return val, "add_float", unit


#
# Functions, may be called by method_func
#
"""
def num_obs_read(dset, field, idx_values):
    dset_read = data.Dataset(dset.rundate, dset.vars['tech'], stage='read', dataset_name=dset.dataset_name,
                             dataset_id='last')
    return _num_obs(dset_read, idx_values)


def num_obs_calculate(dset, field, idx_values):
    dset_calculate = data.Dataset(dset.rundate, dset.vars['tech'], stage='calculate', dataset_name=dset.dataset_name,
                                  dataset_id='last')
    return _num_obs(dset_calculate, idx_values)


def num_obs_estimate(dset, field, idx_values):
    return _num_obs(dset, idx_values)


def rms_estimate(dset, field, idx_values):
    return _rms(dset, idx_values, 'residual')


def rms_calculate(dset, field, idx_values):
    dset_calculate = data.Dataset(dset.rundate, dset.vars['tech'], stage='calculate', dataset_name=dset.dataset_name,
                                  dataset_id='last')
    return _rms(dset_calculate, idx_values, 'residual')


def rms_sisre(dset, field, idx_values):
    return _rms(dset, idx_values, 'sisre')


def rms_sisre_orb(dset, field, idx_values):
    return _rms(dset, idx_values, 'sisre_orb')
"""


def rms(dset, field, idx_values):
    rms = np.array([dset.rms(field, idx=i) for i in _filter_each(dset, idx_values)])
    return rms, "add_float", dset.unit(field)


def num(dset, field, idx_values):
    num_obs = np.array([sum(i) for i in _filter_each(dset, idx_values)])
    return num_obs, "add_float", "count"


def num_clock_breaks(dset, field, idx_values):
    """Number of clock breaks in a dataset
    """
    num_obs = len(idx_values[list(idx_values.keys()).pop()])
    stations = np.array(idx_values.get("station", ("",) * num_obs))

    clock_break_stations = [s for _, _, s in dset.get_events("clock_break")]
    num_clock_breaks = np.array([clock_break_stations.count(s) for s in stations])
    num_clock_breaks[stations == "all"] = len(clock_break_stations)

    return num_clock_breaks, "add_float", "count"


#
# Helper functions
#
def _filter_each(dset, idx_values):
    fields = list(idx_values.keys())
    values = list(zip(*idx_values.values()))

    for value in values:
        filter = {f: v for f, v in zip(fields, value) if v != "all"}
        yield dset.filter(**filter)


def num_obs(dset, idx_values):
    num_obs = np.array([sum(i) for i in _filter_each(dset, idx_values)])
    return num_obs, "add_float", "count"
