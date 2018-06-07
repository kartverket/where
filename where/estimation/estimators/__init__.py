"""Framework for running estimators

Description:
------------

Each estimator should be defined in a separate .py-file. The function inside the .py-file that
should be called need to be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def estimate_least_square(dset):
        ...

The decorated function will be called with a single parameter, ``dset`` which contains a
:class:`~where.data.dataset.Dataset` with data that can be used in the estimation. If the estimator needs the partial
derivatives it should obtain them by itself calling the :func:`where.estimation.partial_vectors` function::

    from where import estimation
    partial_vectors = estimation.partial_vectors(dset, 'estimate', 'estimate_stochastic')




"""

# Where imports
from where.lib import config
from where.lib import plugins


def call(config_key, dset, partial_vectors, obs_noise):
    """Call an estimator

    Args:
        config_key (String):     Config key specifying the name of the estimator.
        dset (Dataset):          Model run data.
        partial_vectors (Dict):  Names and values of the partial derivatives for each partial config key.
        obs_noise (Array):       Observation noise, numpy array with one float value for each observation.
    """
    estimator_name = config.tech[config_key].str
    if estimator_name:
        plugins.call_one(
            package_name=__name__,
            plugin_name=estimator_name,
            dset=dset,
            partial_vectors=partial_vectors,
            obs_noise=obs_noise,
        )


def partial_config_keys(estimator_config_key):
    """Find which partials that a given estimator requires.

    Finds the partial config keys by calling the function registered as 'partial_config_keys' on the given estimator.

    Args:
        estimator_config_key (String):   Config key specifying the name of the estimator.

    Returns:
        Tuple: Strings with names of config keys listing which partial models to run.
    """
    estimator_name = config.tech[estimator_config_key].str
    return plugins.call_one(package_name=__name__, plugin_name=estimator_name, part="partial_config_keys")
