"""Framework for calculating partial derivatives

Description:
------------

Each type of partial derivative should be defined in a separate .py-file. The function inside the .py-file that
should be called need to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midgard.dev import plugins

    @plugins.register
    def calculate_partial_derivative(dset):
        ...

The decorated function will be called with a single parameter, ``dset`` which contains a
:class:`~where.data.dataset.Dataset` with data that can be used when calculating the partial derivatives.

"""
# Standard imports
import datetime

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where.data.time import Time
from where.estimation import estimators
from where.lib import config
from where.lib import log


def partial_vectors(dset, estimator_config_key):
    """Call all partials specified in the configuration and set up the corresponding state vector

    The list of partials to calculate is taken from the config file of the given technique. Each partial calculator is
    passed a :class:`~where.data.dataset.Dataset` with data for the modelrun and should return a tuple with the partial
    vectors and their names.

    Args:
        dset (Dataset):                 A Dataset containing model run data.
        estimator_config_key (String):  Key in config file with the name of the estimator.

    Returns:
        Dict: List of names of the partial derivatives for each partial config key.
    """
    partial_vectors = dict()
    prefix = dset.vars["pipeline"]

    # Delete values from previous iterations
    if "partial" in dset.fields:
        del dset.partial

    for config_key in estimators.partial_config_keys(estimator_config_key):
        partial_vectors[config_key] = list()
        partials = config.tech[config_key].list
        partial_data = plugins.call_all(package_name=__name__, plugins=partials, prefix=prefix, dset=dset)

        for param, (data, names, data_unit) in partial_data.items():
            param_unit_cfg = config.tech[param].unit
            if not param_unit_cfg.str:
                log.fatal(f"No unit given for parameter {param!r} in {param_unit_cfg.source}")

            display_unit = config.tech[param].display_unit.str
            display_unit = param_unit_cfg.str if not display_unit else display_unit
            partial_unit_str = f"{dset.unit('calc')[0]} / ({param_unit_cfg.str})"
            partial_unit = str(Unit(partial_unit_str).u)
            factor = Unit(data_unit, partial_unit)
            data, names = pwlo(config_key, param, dset, data, names)
            for values, name in zip(data.T, names):
                partial_name = f"{param}-{name}" if name else f"{param}"
                partial_vectors[config_key].append(partial_name)

                field_name = f"partial.{partial_name}"
                dset.add_float(field_name, val=values * factor, unit=partial_unit, write_level="operational")
                dset.meta.add(partial_name, display_unit, section="display_units")

    return partial_vectors

def pwlo(config_key, param, dset, data, names):
    """ Set up constant parameters as piecewise local offsets
    
    Estimation epochs and time interval between epochs is defined by the configuration. If no configuration is 
    defined the parameter will be estimated as one constant value the whole timespan of the data.
    
    References:
    Kamil Teke 2011: Sub-daily parameter estimation in {VLBI} data analysis (cite:`teke2011`)
    
    Args:
        config_key (String):    What type of parameter this is. Defined by the configuration file.
        param (String):         Name of the parameter
        dset (Dataset):         Dataset containing model run data.
        data (Array):           Partial derivatives for the given parameter. Size: num_obs x len(names)
        names (List):           Column names for the partial derivatives.
    
    Returns:
        data (Array):           New partial derivatives for the given parameter with pwlo. Size: num_obs x len(new_names) 
        names (List):           New column names for the partial derivatives. New Size: len(old_names)*num_epochs
    """
    if config_key == "estimate_stochastic":
        # Do nothing. pwlo only supported for constant parameters
        return data, names
    
    epoch = config.tech[param].epoch.str
    interval = config.tech[param].knot_interval.str
    if not (epoch and interval):
        # Do not change anything: Estimate one parameter for the whole dataspan with middle of dataspan as epoch
        return data, names
    
    
    partials = []
    column_names = []
    
    for i, n in enumerate(names):
        # Find estimation epochs and column names
        d = data[:, i]
        t_epoch = datetime.datetime.strptime(f"{dset.vars['rundate']} {epoch}", "%Y-%m-%d %H:%M:%S")
        delta = datetime.timedelta(seconds=int(interval))
        t_first = min(dset.time.utc.datetime)
        t_last = max(dset.time.utc.datetime)
        column_name = []
        t = []
                 
        while t_epoch < t_last:
            t_epoch = t_epoch + delta
            if t_epoch >= t_first:
                t.append(t_epoch)
                column_name.append(f"{n}:{t_epoch:%Y-%m-%dT%H:%M:%S}")
        
        t_epoch = t[0] - delta
        t.insert(0, t_epoch)
        column_name.insert(0, f"{n}:{t_epoch:%Y-%m-%dT%H:%M:%S}")
        
        t = Time(t, scale="utc", fmt="datetime").mjd
        
        # TODO: Make sense in timeseries?

        d_dx = np.zeros((dset.num_obs, len(column_name)))  
    
        for j in range(len(column_name) - 1):
            dx_d1 = np.zeros(dset.num_obs)
            dx_d2 = np.zeros(dset.num_obs)
            idx = np.logical_and(dset.time.mjd > t[j], dset.time.mjd < t[j+1])
            
            dx_d1[idx] = 1 - (dset.time.mjd[idx] - t[j])/(t[j+1] - t[j])
            dx_d2[idx] = (dset.time.mjd[idx] - t[j])/(t[j+1] - t[j])
            
            d_dx[idx, j] = dx_d1[idx]
            d_dx[idx, j + 1]= dx_d2[idx]
                     
        partials.append(d[:, None] * d_dx)
        column_names.append(column_name)
    
    return np.concatenate(partials, axis=1), list(np.concatenate(column_names))

    
    