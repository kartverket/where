"""Framework for calculating models

Description:
------------

Each model should be defined in a separate .py-file inside one of the model directories, delay, site, orbit and
satellite. The function inside the .py-file that should be called need to be decorated with the
:func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def vlbi_vacuum_delay(dset):
        ...

See the corresponding :file:`__init__.py`-files for more information.


See also:
---------

* :mod:`where.models.delay`
* :mod:`where.models.site`

"""
# External library imports

# Where imports
from where.lib import log
from where.models import delay
from where.models import site

# Make documentation functions available
from where.models.delay import doc as doc_delay  # noqa
from where.models.site import doc as doc_site  # noqa

# Make orbit functions optionally available (does not require compilation if orbits are not used)
try:
    from where.models.orbit._orbit import calculate as calculate_orbit, update_orbit  # noqa
except ImportError:
    from midgard.dev import optional

    calculate_orbit = optional.SimpleMock(
        "where.models.orbit._orbit.calculate", error_msg="Try running setup_cython.py to compile the orbit modules"
    )
    update_orbit = optional.SimpleMock(
        "where.models.orbit._orbit.update_orbit", error_msg="Try running setup_cython.py to compile the orbit modules"
    )


def calculate_delay(config_key, dset_in, dset_out=None, shape=(), write_levels=None):
    """Call delay models and store output in dataset

    Args:
        config_key (String):  Key in config with list of models.
        dset_in (Dataset):    Dataset to read data from.
        dset_out (Dataset):   Dataset to store data to.
        shape (Tuple of int): Shape of output.
    """
    _calculate_model(delay.calculate, config_key, dset_in, dset_out, shape, write_levels)


def calculate_site(config_key, dset_in, dset_out=None, shape=(3,), write_levels=None):
    """Call site models and store output in dataset

    Args:
        config_key (String):  Key in config with list of models.
        dset_in (Dataset):    Dataset to read data from.
        dset_out (Dataset):   Dataset to store data to.
        shape (Tuple of int): Shape of output.
    """
    _calculate_model(site.calculate, config_key, dset_in, dset_out, shape, write_levels)


def _calculate_model(calculate_func, config_key, dset_in, dset_out, shape, write_levels=None):
    """Call models and store output in dataset

    If the model output is empty, we still create a dummy field in the table only containing zeros. This is done to
    assert that the table will always exist after doing a `models.calculate...`-call.

    Args:
        calculate_func (Function):  The function that calls models.
        config_key (String):        Key in config with list of models, also table the model output is stored in.
        dset_in (Dataset):          Dataset to read data from.
        dset_out (Dataset):         Dataset to store data to.
        shape (Tuple of int):       Shape of output.
    """
    dset_out = dset_in if dset_out is None else dset_out
    write_levels = dict() if write_levels is None else write_levels

    model_output = calculate_func(config_key, dset_in)
    if not model_output:
        model_name = "{}_zeros".format(config_key)
        if model_name not in dset_out.fields:
            dset_out.add_float(
                model_name,
                table=config_key,
                shape=shape,
                unit="meter",
                write_level=write_levels.get(model_name, "analysis"),
            )

    for model_name, values in sorted(model_output.items()):
        if model_name in dset_out.fields:
            # TODO: Add possibilty to change table to 'calc_models' from a 'float_data' table 
            #       (dset_out._fields[model_name] = config_key)
            dset_out[model_name][:] = values
        else:
            dset_out.add_float(
                model_name,
                table=config_key,
                shape=shape,
                val=values,
                unit="meter",
                write_level=write_levels.get(model_name, "analysis"),
            )

        log.info(f"Average correction = {dset_out.rms(model_name):14.5f} in {model_name} model")
