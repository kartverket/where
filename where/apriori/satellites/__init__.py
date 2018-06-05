"""Framework for satellites

Description:
------------

To add a new satellite, simply create a new .py-file which defines a class inheriting from satellites.Satellite. The
class needs to be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.apriori.satellites import satellite
    from where.lib import plugins

    @plugins.register
    class Satellite_1(satellite.Satellite):
        ...

To use a satellite, you will typically use the :func:`get_satellite`-function defined below

    from where import apriori
    my_new_satellite = apriori.get_satellite('satellite_1', ...)

The name used in `get` to call the satellite is the lower case name of the class containing the satellite info.




$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""

# Where imports
from where.lib import cache
from where.lib import log
from where.lib import plugins

# Do not support * imports
__all__ = []


@cache.function
def satellites():
    """Available satellites and where they are defined

    Returns:
        Dict: For each satellite, a tuple describing which satellite plug-in defines the satellite.
    """
    return {p.lower(): (s, p) for s in plugins.list_all(__name__) for p in plugins.list_parts(__name__, s)}


@cache.function
def names():
    """List the names of the available satellites

    Returns:
        List: List of strings with the names of the available satellites
    """
    return sorted(satellites().keys())


@cache.function
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
        log.fatal("Unknown satellite '{}'. Defined satellites are {}.", satellite_name, ", ".join(names()))

    return plugins.call_one(package_name=__name__, plugin_name=plugin, part=part, **kwargs)


@cache.function
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
