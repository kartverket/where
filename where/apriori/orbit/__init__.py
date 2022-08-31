"""Framework for reading and working with initial (precalculated) orbits

Description:
------------

Each type of initial orbit should be defined as a class in a separate .py-file. The class inside the .py-file that
represents the initial orbit need to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midgard.dev import plugins

    @plugins.register
    class BroadcastOrbit(orbit.InitialOrbit):
        ...

The initial orbit class will be responsible for parsing the initial orbit, as well as doing necessary calculations and
interpolations. It must inherit from apriori.orbit.InitialOrbit or a subclass.

To obtain a initial orbit, call the :func:`get_orbit`-function. Which initial orbit to use can be specified explicitly. 
If it is not, the initial orbit is looked up in the config files. See :func:`get_orbit` for more information.
"""
# Standard library imports
from functools import lru_cache

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log

# Make available for convenience
from where.apriori.orbit._orbit import AprioriOrbit  # noqa


@lru_cache()
def get_orbit(apriori_orbit=None, **kwargs):
    """Get an apriori orbit

    The specification of the apriori orbit is matched with a filename in this orbit-directory. If it is not passed in
    as an argument, the apriori orbit to use is read from the configuration.

    Args:
        apriori_orbit (String):  Optional specification of which apriori orbit to use (see above).

    Returns:
        AprioriOrbit: Apriori orbit object.
    """
    apriori_orbit = config.tech.get("apriori_orbit", apriori_orbit).str
    if apriori_orbit not in ["has", "broadcast", "precise", "slr"]:
        log.fatal(
            f"Configuration value '{apriori_orbit}' for option 'apriori_orbit' is unknown. It should be either 'broadcast' "
            "and/or 'precise', or 'has' or 'slr'. "
        )

    return plugins.call(package_name=__name__, plugin_name=apriori_orbit, **kwargs)
