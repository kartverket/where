"""Framework for reading apriori data sources

Description:
------------

Each data source should be defined in a separate .py-file. The function inside the .py-file that
should be called need to be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def read_fun_datasource(rundate):
        ...

The decorated function will be called through the :func:`get`. Parameters to the registered function should be passed
as named keyword arguments.




"""
# Where imports (more imports are done locally to avoid circular imports)
from where.lib import cache
from where.lib import exceptions
from where.lib import log
from where.lib import plugins

# Make functions in subpackages available as apriori-plugins
from where.apriori.orbit import get_orbit  # noqa
from where.apriori.satellites import get_satellite, get_satellite_vars  # noqa
from where.apriori.trf import get_trf
from where.apriori.crf import get_crf

plugins.register(get_orbit)
plugins.register(get_satellite)
plugins.register(get_trf)
plugins.register(get_crf)


def names():
    """List the names of the available parsers

    Returns:
        List: List of strings with the names of the available parsers
    """
    return plugins.list_all(package_name=__name__)


@cache.function
def get(datasource_name, **kwargs):
    """Read data from the given data source

    Simple data sources that only return data directly from a parser does not need an explicit apriori-file. This is
    handled by looking in the parser-directory if a data source is not found in the apriori directory.

    The import of where.parsers is done locally to avoid circular imports.

    Args:
        datasource_name (String):   Name of apriori data source
        kwargs:                     Input arguments to the data source

    Returns:
        The data from the data source (data type depends on source)
    """
    try:
        return plugins.call_one(package_name=__name__, plugin_name=datasource_name, **kwargs)
    except exceptions.UnknownPluginError as apriori_err:
        from where import parsers

        try:
            data = parsers.parse_key(file_key=datasource_name, **kwargs).as_dict()
            log.dev(f"Called parsers.parse_key({datasource_name}) in apriori.get()")
            return data
        except (AttributeError) as att:
            try:
                data = parsers.parse(datasource_name, **kwargs)
                log.dev(f"Called parsers.parse({datasource_name}) in apriori.get()")
                return data
            except exceptions.UnknownPluginError:
                raise apriori_err from None
