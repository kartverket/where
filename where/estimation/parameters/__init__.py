"""Framework for calculating partial derivatives

Description:
------------

Each type of partial derivative should be defined in a separate .py-file. The function inside the .py-file that
should be called need to be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def calculate_partial_derivative(dset):
        ...

The decorated function will be called with a single parameter, ``dset`` which contains a
:class:`~where.data.dataset.Dataset` with data that can be used when calculating the partial derivatives.




"""

# Where imports
from where.lib import config
from where.estimation import estimators
from where.lib import log
from where.lib import plugins
from where.lib.unit import unit


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
    prefix = config.analysis.get("analysis", default="").str

    for config_key in estimators.partial_config_keys(estimator_config_key):
        partial_vectors[config_key] = list()
        partial_data = plugins.call_all(package_name=__name__, config_key=config_key, prefix=prefix, dset=dset)

        for param, (data, names, data_unit) in partial_data.items():
            param_unit_cfg = config.tech[param].unit
            log.assert_true(param_unit_cfg.str, "No unit given for parameter '{}' in {}", param, param_unit_cfg.source)

            display_unit = config.tech[param].display_unit.str
            display_unit = param_unit_cfg.str if not display_unit else display_unit
            partial_unit = str(unit("{} / ({})".format(dset.unit("calc"), param_unit_cfg.str)).u)
            factor = unit(data_unit, partial_unit)
            for values, name in zip(data.T, names):
                partial_name = "{}-{}".format(param, name)
                dset.add_float("partial_" + partial_name, table="partial", val=values * factor, unit=partial_unit)
                dset.add_to_meta("display_units", partial_name, display_unit)
                partial_vectors[config_key].append(partial_name)

    return partial_vectors
