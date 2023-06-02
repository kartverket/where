"""Write key indicators from model run to timeseries dataset

Description:
------------

We store some indicators from a daily analysis to a common dataset with a dummy-date of January 1st 1970. This is
called the timeseries dataset and can be used to look at results across different datasets.

"""

# Standard library imports
from datetime import date
import itertools
import re
import sys

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.data import dataset3 as dataset
from where.lib import log


WRITER = __name__.split(".")[-1]


@plugins.register
def add_to_full_timeseries(dset):
    """Write some key variables to the full timeseries

    Args:
        dset:  Dataset, data for a model run.
    """
    log.info(f"Updating timeseries dataset")

    dset_session = dataset.Dataset()

    # Add data to dset_session
    idx_fields = config.tech[WRITER].index.list
    field_values = []
    for f in idx_fields:
        dset_value = set(dset.unique(f)) # field with data in dataset
        meta_value = set(dset.meta[f].keys()) # field without data in dataset, but is defined in meta
        values = dset_value.union(meta_value)
        field_values.append(["all"] + sorted(list(values)))
    #field_values = [set("all") + set(dset.unique(f)) + set(dset.meta[f].keys()) for f in idx_fields]
    idx_values = dict(zip(idx_fields, zip(*itertools.product(*field_values))))
    # TODO: Remove combinations where filter leaves 0 observations

    num_obs = len(idx_values[idx_fields[0]])  # Length of any (in this case the first) field
    mean_epoch = dset.time.mean.utc
    rundate_str = dset.analysis["rundate"].strftime(config.FMT_date)
    session_code = dset.vars.get("session_code", "")
    status = dset.meta.get("analysis_status", "unchecked")
    session_type = dset.meta.get("input", dict()).get("where_session_type", "")

    dset_session.num_obs = num_obs
    dset_session.add_time("time", val=[mean_epoch] * num_obs, scale=mean_epoch.scale, fmt=mean_epoch.fmt)
    dset_session.add_text("rundate", val=[rundate_str] * num_obs)
    dset_session.add_text("status", val=[status] * num_obs)
    dset_session.add_text("session_type", val=[session_type] * num_obs)
    dset_session.add_text("session_code", val=[session_code] * num_obs)

    for field, value in idx_values.items():
        dset_session.add_text(field, val=value)

    default_dset_str = f"{dset.vars['stage']}/{dset.vars['label']}"
    dsets = {default_dset_str: dset}
    for method, cfg_entry in config.tech[WRITER].items():
        try:
            method_func = getattr(sys.modules[__name__], f"method_{method}")
        except AttributeError:
            log.warn(f"Method {method!r} is unknown")
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
                dsets[dset_str] = dataset.Dataset.read(
                    rundate=dset.analysis["rundate"],
                    pipeline=dset.vars["pipeline"],
                    stage=stage,
                    session_code=dset.vars["session_code"],
                    label=dset_id,
                    id=dset.analysis["id"],
                )

            val, adder, unit = method_func(dsets[dset_str], field_in, idx_values, func)
            if adder:
                add_func = getattr(dset_session, adder)
                add_func(field_out, val=val) if unit is None else add_func(field_out, val=val, unit=unit)
                

    # hack to get solved neq data into the time series:
    # TODO: unhack this :P Add as a method_neq instead?
    if "normal equation" in dset.meta:
        _add_solved_neq_fields(dset, dset_session, idx_values)

    # Read timeseries dataset and extend it with session dataset
    dset_id = config.tech.timeseries.dataset_id.str.format(**dset.vars)
    try:
        # Read existing dataset
        dset_ts = dataset.Dataset.read(
            rundate=date(1970, 1, 1),
            pipeline=dset.vars["pipeline"],
            stage="timeseries",
            label=dset_id,
            session_code="",
            use_options=False,
            id=dset.analysis["id"],
        )
    except OSError:
        # Start new timeseries dataset
        dset_ts = dataset.Dataset(
            rundate=date(1970, 1, 1),
            pipeline=dset.vars["pipeline"],
            stage="timeseries",
            label=dset_id,
            session_code="",
            use_options=False,
            id=dset.analysis["id"],
        )

    if dset_ts.num_obs > 0:
        # Filter timeseries dataset to remove any previous data for this rundate and session
        keep_idx = np.logical_not(dset_ts.filter(rundate=rundate_str, session_code=session_code))
        dset_ts.subset(keep_idx)

    # Extend timeseries dataset with dset_session and write to disk
    dset_ts.extend(dset_session)
    dset_ts.write()


def _add_solved_neq_fields(dset, dset_session, idx_values):
    names = dset.meta["normal equation"]["names"]
    x = np.array(dset.meta["normal equation"]["solution"])
    Q = np.array(dset.meta["normal equation"]["covariance"])
    units = np.array(dset.meta["normal equation"]["unit"])

    # This code is terrible. TODO: Rewrite
    idx_names = idx_values[list(idx_values.keys()).pop()]
    num_obs = len(idx_names)
    fields = np.unique([n.split("-")[0] for n in names])
    for field in fields:
        if "vlbi_src_dir" in field or "vlbi_baseline" in field:
            # TODO
            continue
        params = [f for f in names if f.startswith(field + "-")]
        param_units = [u for f, u in zip(names, units) if f.startswith(field + "-")]
        name, name2 = _parse_params(params, idx_names)
        epochs = [n.find(":") > -1 for n in name2]
        if any(epochs):
            log.warn(f"NEQ parameter ({field}) with estimation epochs different from mid-session cannot be added to timeseries")
            continue
    
        dim = len(np.unique(name2)) if name2 else 1
        shape = (num_obs, dim) if dim > 1 else num_obs
        shape_cov = (num_obs, dim, dim) if dim > 1 else num_obs
        val = np.full(shape, np.nan, dtype=float)
        val_cov = np.full(shape_cov, np.nan, dtype=float)
        idx = np.array([any([n == p for p in name]) for n in idx_names], dtype=bool)
        mean = np.array([state for state, param in zip(x, names) if param in params]).reshape(-1, dim)

        if idx.any():
            val[idx] = mean.reshape(val[idx].shape)
            for i, n in enumerate(idx_names):
                if not idx[i]:
                    continue
                idx_param = np.array(
                    [True if param in params and n == param.split("-", maxsplit=1)[1].rsplit("_", maxsplit=1)[0] else False for param in names], dtype=bool
                )
                param_cov = Q[idx_param][:, idx_param]
                val_cov[i] = param_cov
        else:
            val[0] = mean
            idx_param = np.array([True if param in params else False for param in names], dtype=bool)
            mean_cov = Q[idx_param][:, idx_param]
            val_cov[0] = mean_cov

        # all fields of the same type share the same unit, only the first entry of param_units is needed
        dset_session.add_float("neq_" + field, val=val, unit=param_units[0])
        dset_session.add_float("neq_" + field + "_cov_", val=val_cov, unit=param_units[0])

    # Add solution helmert parameters to timeseries
    sections = "helmert", "neq_helmert"
    fields = ["T_X", "T_Y", "T_Z", "scale", "alpha", "beta", "gamma"]
    for section in sections:
        for field in fields:
            try:
                value, unit = dset.meta[section][field]
                data = np.full(num_obs, np.nan, dtype=float)
                data[0] = value
                dset_session.add_float(f"{section}_{field}", val=data, unit=unit)
            except KeyError:
                pass


#
# Methods, corresponds to keys in config.tech.timeseries
#
# TODO: Use config metadata instead of dummy functions below
def method_index(dset, field, idx_values, func):
    return None, None, None


def method_dataset_name(dset, field, idx_values, func):
    return None, None, None


def method_dataset_id(dset, field, idx_values, func):
    return None, None, None


def method_func(dset, field, idx_values, func):
    calculate_func = getattr(sys.modules[__name__], func)
    return calculate_func(dset, field, idx_values)


def method_statistics(dset, field, idx_values, func):
    num_obs = len(idx_values[list(idx_values.keys()).pop()])
    val = np.full(num_obs, np.nan, dtype=float)
    val[:] = dset.meta["statistics"][field]

    return val, "add_float", None  # Todo: Can we add a unit here somehow?

def method_meta_index(dset, field, idx_values, func):
    
    key = list(idx_values.keys()).pop()
    names = idx_values[key]
    num_obs = len(names)
    val = np.full(num_obs, np.nan, dtype=float)
    for i, n in enumerate(names):
        if n == "all":
            continue
        val[i] = dset.meta[key][n].get(field, np.nan)
        
    return val, "add_float", None

def method_meta(dset, field, idx_values, func):
    """Save simple meta information in timeseries
    
    Fields must be saved in meta in the following way:
        dset.meta[field] = {'value': 146.41794074611778, '__unit__': 'Megameter**3'}

    "value" is assumed to be a float
    """
    num_obs = len(idx_values[list(idx_values.keys()).pop()])
    val = np.full(num_obs, np.nan, dtype=float)
    val[:] = dset.meta[field]["value"]
    unit = dset.meta[field].get("__unit__")
    
    return val, "add_float", unit


def method_text(dset, field, idx_values, func):
    if field not in dset.fields:  # Ignore fields not in dataset
        return None, None, None

    text = np.array(["/".join(dset.unique(field, **f)) for f in _filter_each(idx_values)])
    return text, "add_text", None


def method_state(dset, field, idx_values, func):
    if "vlbi_src_dir" in field:
        # TODO
        return None, None, None
    idx_names = idx_values[list(idx_values.keys()).pop()]
    num_obs = len(idx_names)
    params = [f for f in dset.fields if f.startswith("state." + field + "-") and not f.endswith("_sigma")]
    name, name2 = _parse_params(params, idx_names)
    epochs = [n.find(":") > -1 for n in name2]
    if any(epochs):
        log.warn(f"Parameter ({field}) with estimation epochs different from mid-session cannot be added to timeseries")
        return None, None, None
    dim = len(np.unique(name2)) if name2 else 1
    val = np.full((num_obs, dim), np.nan, dtype=float)
    idx = np.array([any([n == p for p in name]) for n in idx_names], dtype=bool)
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

def _parse_params(params, idx_names):
    "Split the parameter name of the subparameter into two parts"
    subparams = [p.split("-", maxsplit=1)[-1] for p in params]
    name = []
    name2 = []

    # Parameters connected to the idx_values, e.g. station params
    for idx_name in idx_names:
        for param in subparams:
            idx = param.find(idx_name) # Assumption: If found idx is always 0.
            if idx >= 0:
                if len(idx_name) < len(param) and param[idx+len(idx_name)] == "_":
                    # Style of param: {station}_{subparam}
                    name.append(param[idx:idx+len(idx_name)])
                    if len(idx_name) < len(param):
                        name2.append(param[idx+len(idx_name):])
                elif len(idx_name) == len(param):
                    # Style of param: {station}
                    name.append(param[idx:idx+len(idx_name)])
    # Parameters not connected the idx_values, e.g. global params
    if not name:
        name = [p.split("_", maxsplit=1)[0] for p in subparams]
        try:
            name2 = [p.split("_", maxsplit=1)[-1] for p in subparams]
        except IndexError:
            pass

    return name, name2

#
# Functions
#


def rms(dset, field, idx_values):
    if field not in dset.fields:  # Ignore fields not in dataset
        return None, None, None

    rms = np.array([dset.rms(field, **f) for f in _filter_each(idx_values)])
    return rms, "add_float", dset.unit(field)


def nanrms(dset, field, idx_values):
    if field not in dset.fields:  # Ignore fields not in dataset
        return None, None, None

    not_nan = ~np.isnan(dset[field])
    rms = np.array([dset.rms(field, idx=not_nan, **f) for f in _filter_each(idx_values)])
    return rms, "add_float", dset.unit(field)


def num(dset, _, idx_values):
    try:
        num_obs = np.array([dset.num(**f) for f in _filter_each(idx_values)])
        return num_obs, "add_float", "count"
    except AttributeError:
        return None, None, None


def count(dset, field, idx_values):
    try:
        counts = np.array([dset.count(field, **f) for f in _filter_each(idx_values)])
        return counts, "add_float", "count"
    except AttributeError:
        return None, None, None


def num_clock_breaks(dset, field, idx_values):
    """Number of clock breaks in a dataset
    """
    num_obs = len(idx_values[list(idx_values.keys()).pop()])
    stations = np.array(idx_values.get("station", ("",) * num_obs))

    clock_break_stations = [s for _, _, s in dset.meta.get_events("clock_break")]
    num_clock_breaks = np.array([clock_break_stations.count(s) for s in stations])
    num_clock_breaks[stations == "all"] = len(clock_break_stations)

    return num_clock_breaks, "add_float", "count"


#
# Helper functions
#


def _filter_each(idx_values):
    fields = list(idx_values.keys())
    values = list(zip(*idx_values.values()))

    for value in values:
        yield {f: v for f, v in zip(fields, value) if v != "all"}
