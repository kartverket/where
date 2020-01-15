"""Framework for integrators

Description:
------------

Each integrator should be defined in a separate .py-file. The function inside the .py-file that should be called
needs to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midard.dev import plugins

    @plugins.register
    def cowell(dset):
        ...

"""

import importlib

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import log

# Do not support * imports
__all__ = []


def names():
    """List names of integrators for a given session

    The list of names will be sorted according to the sort_value of the integrator plugins (use `register_ordered` to
    define a sort value for an integrator).

    Returns:
        List:   Sorted list of names of integrators.
    """
    return plugins.names(package_name=__name__)


def call(integrator, **kwargs):
    """Call one integrator

    Args:
        integrator (String):   The integrator to call.
        **kwargs:              Arguments to be passed to function.

    Return:
        Output of integrator
    """
    try:
        integrator_module = importlib.import_module(f"{__name__}.{integrator}")
        return integrator_module.integrate(**kwargs)
    except ImportError:
        log.fatal(f"Did not find compiled Cython module '{__name__}.{integrator}'")
