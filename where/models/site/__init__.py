"""Framework for calculating site displacement models

Description:
------------

Each site displacement model should be defined in a separate .py-file. The function inside the .py-file that should be
called need to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midgard.dev import plugins

    @plugins.register
    def solid_tides(dset):
        ...

The decorated function will be called with a single parameter, ``dset`` which contains a
:class:`~where.data.dataset.Dataset` with data that can be used when calculating the site displacement.


"""
# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log
from where.data import position


def calculate(config_key, dset):
    prefix = dset.vars["pipeline"]
    return plugins.call_all(package_name=__name__, plugins=config.tech[config_key].list, prefix=prefix, dset=dset)


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
        for multiplier, pos_delta in zip(dset_out.for_each_suffix("station"), values):
            field_name = f"{config_key}.{model_name}{dset_out.default_field_suffix}"
            if field_name in dset_out.fields:
                dset_out[field_name][:] = pos_delta
            else:
                if position.is_position_delta(pos_delta):
                    dset_out.add_position_delta(
                        field_name,
                        pos_delta,
                        multiplier=multiplier,
                        write_level=write_levels.get(model_name, "analysis"),
                    )
                elif position.is_posvel_delta(pos_delta):
                    dset_out.add_posvel_delta(
                        field_name,
                        pos_delta,
                        multiplier=multiplier,
                        write_level=write_levels.get(model_name, "analysis"),
                    )
        log.info(f"Average correction = {dset_out.rms(f'{config_key}.{model_name}'):14.5f} in {model_name} model")


def add(config_key, dset):
    """Sum all site model corrections based on caculated site model field results in GCRS

    Args:
        config_key (String):  Key in config with list of models.
        dset (Dataset):       A dataset containing the data.

    Returns:
        List with PositionDelta object array with total sum of site corrections given for each station entry
    """
    delta_pos = list()

    for _ in dset.for_each_suffix("station"):
        if position.is_position(dset.site_pos):
            sta_delta = position.PositionDelta(np.zeros((dset.num_obs, 3)), system="gcrs", ref_pos=dset.site_pos)
        elif position.is_posvel(dset.site_pos):
            sta_delta = position.PosVelDelta(np.zeros((dset.num_obs, 6)), system="gcrs", ref_pos=dset.site_pos)

        if config_key not in dset.fields:
            return list(sta_delta)

        pos_fields = [f for f in dset[config_key]._fields if f.endswith(dset.default_field_suffix)]

        for field in pos_fields:
            sta_delta += dset[config_key][field]
        delta_pos.append(sta_delta)

    return delta_pos


def calculate_site(config_key, dset_in, dset_out=None, write_levels=None):
    """Call site models and store output in dataset

    Args:
        config_key (String):  Key in config with list of models.
        dset_in (Dataset):    Dataset to read data from.
        dset_out (Dataset):   Dataset to store data to.
    """
    _calculate_model(calculate, config_key, dset_in, dset_out, write_levels)
