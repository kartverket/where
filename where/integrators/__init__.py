"""Framework for integrators

Description:
------------

Each integrator should be defined in a separate .py-file. The function inside the .py-file that should be called
needs to be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def cowell(dset):
        ...



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

import importlib

# Where imports
from where.lib import plugins, log


# Do not support * imports
__all__ = []


def names():
    """List names of integrators for a given session

    The list of names will be sorted according to the sort_value of the integrator plugins (use `register_ordered` to
    define a sort value for an integrator).

    Returns:
        List:   Sorted list of names of integrators.
    """
    return plugins.list_all(package_name=__name__)


def call(integrator, **kwargs):
    """Call one integrator

    Args:
        integrator (String):   The integrator to call.
        **kwargs:              Arguments to be passed to function.

    Return:
        Output of integrator
    """
    try:
        integrator_name = importlib.import_module("{}.{}".format(__name__, integrator))
        # integrator_function = __name__.rsplit('.', maxsplit=1)[1]
        return integrator_name.integrate(**kwargs)
    except ImportError:
        log.fatal("Did not find compiled Cython module '{}.{}'", package, integrator_name)
