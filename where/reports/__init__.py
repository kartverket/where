"""Framework for writing reports from Where model runs

Description:
------------

Each report should be defined in a separate .py-file. The function inside the .py-file that should be called need to be
decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def write_cool_report(dset):
        ...

The decorated function will be called with a single parameter, ``dset`` which contains a
:class:`~where.data.dataset.Dataset` with the data that should be output.




"""

# Where imports
from where.lib import config
from where.lib import plugins

# Make report available as a module on the report package, so we can use reports.report.
from where.reports import report  # noqa


def names():
    """List the names of the available reports specified in the configuration

    The available reports (modules in the reports directory) are compared to the list specified in the ``report``-field
    of the configuration file. Only reports that appears both places are returned.

    Returns:
        List: List of strings with the names of the available reports.

    """
    return plugins.list_all(package_name=__name__, config_key="report")


def write(rundate, tech):
    """Call all reports specified in the configuration

    The list of reports to use is taken from the config file of the given technique. Each report is passed a
    :class:`~where.data.dataset.Dataset` with data for the modelrun and should write the relevant parts of the data to
    file.

    Args:
        data:     A Dataset containing model run data.

    """
    prefix = config.analysis.get("analysis", default="").str
    plugins.call_all(package_name=__name__, config_key="report", prefix=prefix, rundate=rundate, tech=tech)
