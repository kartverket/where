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

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where.lib import config
from where.estimation import estimators
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
            for values, name in zip(data.T, names):
                partial_name = f"{param}-{name}" if name else f"{param}"
                partial_vectors[config_key].append(partial_name)

                field_name = f"partial.{partial_name}"
                dset.add_float(field_name, val=values * factor, unit=partial_unit, write_level="operational")
                dset.meta.add(partial_name, display_unit, section="display_units")

    return partial_vectors
