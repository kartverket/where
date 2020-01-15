"""Framework for constructing a Dataset from observation files and apriori sources

Description:
------------

Each pipeline should be defined separate .py-file. The function inside the .py-file that
should be called need to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midgard_dev import plugins

    @plugins.register
    def write_to_dataset(dset, rundate=None, session=None, **obsargs):
        ...

The decorated function will be called through the :func:`get`, with rundate, pipeline and session specificed as input
arguments. Additional parameters to the registered function should be passed as named keyword arguments.
"""

# Midgard imports
from midgard.dev import plugins


def get(dset, **obs_args):
    """Construct a Dataset for the pipeline based on observations

    Args:
        dset:  A Dataset that will be filled with observations and necessary fields
    """
    rundate = dset.analysis["rundate"]
    pipeline = dset.vars["pipeline"]

    plugins.call(package_name=__name__, plugin_name=pipeline, dset=dset, rundate=rundate, **obs_args)
