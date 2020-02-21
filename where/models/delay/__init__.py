"""Framework for calculating delay models

Description:
------------

Each delay model should be defined in a separate .py-file. The function inside the .py-file that should be called need
to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midgard.dev import plugins

    @plugins.register
    def vlbi_vacuum_delay(dset):
        ...

The decorated function will be called with a single parameter, ``dset`` which contains a
:class:`~where.data.dataset.Dataset` with data that can be used when calculating the delay.

"""
# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log


def calculate_delay(config_key, dset_in, dset_out=None, write_levels=None):
    """Call delay models and store output in dataset

    Args:
        config_key (String):  Key in config with list of models.
        dset_in (Dataset):    Dataset to read data from.
        dset_out (Dataset):   Dataset to store data to.
    """
    _calculate_model(calculate, config_key, dset_in, dset_out, write_levels)


def _calculate_model(calculate_func, config_key, dset_in, dset_out, write_levels=None):
    """Call models and store output in dataset

    If the model output is empty, we still create a dummy field in the table only containing zeros. This is done to
    assert that the table will always exist after doing a `models.calculate...`-call.

    Args:
        calculate_func (Function):  The function that calls models.
        config_key (String):        Key in config with list of models, also table the model output is stored in.
        dset_in (Dataset):          Dataset to read data from.
        dset_out (Dataset):         Dataset to store data to.
    """
    dset_out = dset_in if dset_out is None else dset_out
    write_levels = dict() if write_levels is None else write_levels

    model_output = calculate_func(config_key, dset_in)

    for model_name, values in sorted(model_output.items()):
        field_name = f"{config_key}.{model_name}"
        if field_name in dset_out.fields:
            dset_out[field_name][:] = values
        else:
            dset_out.add_float(field_name, values, write_level=write_levels.get(model_name, "analysis"), unit="meter")
        log.info(f"Average correction = {dset_out.rms(f'{config_key}.{model_name}'):14.5f} in {model_name} model")


def add(config_key, dset):
    delay_fields = [f for f in dset[config_key]._fields]
    delta_delay = np.zeros(dset.num_obs)
    for field in delay_fields:
        delta_delay += dset[config_key][field]

    return delta_delay


def calculate(config_key, dset):
    prefix = dset.vars["pipeline"]
    return plugins.call_all(package_name=__name__, plugins=config.tech[config_key].list, prefix=prefix, dset=dset)
