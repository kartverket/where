"""Framework for satellites

Description:
------------

To add a new satellite, simply create a new .py-file which defines a class inheriting from satellites.Satellite. The
class needs to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from where.apriori.satellites import satellite
    from midgard.dev import plugins

    @plugins.register
    class Satellite_1(satellite.Satellite):
        ...

To use a satellite, you will typically use the :func:`get_satellite`-function defined below

    from where import apriori
    my_new_satellite = apriori.get_satellite('satellite_1', ...)

The name used in `get` to call the satellite is the lower case name of the class containing the satellite info.
"""
# Standard library imports
from functools import lru_cache

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import log

# Do not support * imports
__all__ = []


@lru_cache()
def satellites():
    """Available satellites and where they are defined

    Returns:
        Dict: For each satellite, a tuple describing which satellite plug-in defines the satellite.
    """
    return {p.lower(): (s, p) for s in plugins.names(__name__) for p in plugins.parts(__name__, s)}


@lru_cache()
def names():
    """List the names of the available satellites

    Returns:
        List: List of strings with the names of the available satellites
    """
    return sorted(satellites().keys())


@lru_cache()
def get_satellite(satellite_name, **kwargs):
    """Get a satellite object by name

    Args:
        satellite_name (String): Name used to look up satellite.
        kwargs (Dict):           Arguments that will be passed to the satellite object.

    Returns:
        A satellite object describing the satellite.
    """
    try:
        plugin, part = satellites()[satellite_name.lower()]
    except KeyError:
        log.fatal(f"Unknown satellite '{satellite_name}'. Defined satellites are {', '.join(names())}")

    return plugins.call(package_name=__name__, plugin_name=plugin, part=part, **kwargs)


@lru_cache()
def get_satellite_vars(satellite_name, **kwargs):
    """Construct a dict of satellite variables

    Args:
        satellite_name (String): Name used to look up satellite.
        kwargs (Dict):           Arguments that will be passed to the satellite object.

    Returns:
        Dict: Dictionary with variables containing basic information about the satellite.
    """
    satellite = get_satellite(satellite_name, **kwargs)

    # Create the dict of satellite variables
    return dict(sat_name=satellite.name, sat_shortname=satellite.short_name, sat_id=str(satellite.cospar_id))
