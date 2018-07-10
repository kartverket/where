"""Framework for writing output in different formats

Description:
------------

Each output format / output destination should be defined in a separate .py-file. The function inside the .py-file that
should be called need to be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def write_as_fancy_format(dset):
        ...

The decorated function will be called with a single parameter, ``dset`` which contains a
:class:`~where.data.dataset.Dataset` with the data that should be output.




"""

# Where imports
from where.lib import config
from where import data
from where.lib import plugins


def names():
    """List the names of the available writers specified in the configuration

    The available writers (modules in the writers directory) are compared to the list specified in the ``output``-field
    of the configuration file. Only writers that appears both places are returned.

    Returns:
        List: List of strings with the names of the available writers.

    """
    return plugins.list_all(package_name=__name__, config_key="output")


def write(default_stage):
    """Call all writers specified in the configuration

    The list of writers to use is taken from the config file of the given technique. Each writer is passed a
    :class:`~where.data.dataset.Dataset` with data for the modelrun and should write the relevant parts of the data to
    file.

    By default the last dataset for the default_stage is sent to the writer, but that is possible to override with the
    following notation:

        output = writer_1                 # Use last dataset of default_stage
        output = writer_1:calculate       # Use last dataset of "calculate" stage
        output = writer_1:calculate/2     # Use dataset 2 of "calculate" stage

    Args:
        default_stage (String):    Name of stage to read dataset from by default.
    """
    dsets = dict()
    prefix = config.analysis.get("analysis", default="").str
    output_list = config.tech.output.list
    writer_and_dset = [f"{o}:{default_stage}".split(":")[:2] for o in output_list]

    rundate = config.analysis.rundate.date
    tech = config.analysis.tech.str
    session = config.analysis.session.str
    for writer, dset_str in writer_and_dset:

        # Read the datasets
        if dset_str not in dsets:
            stage, _, dset_id = dset_str.partition("/")
            dset_id = int(dset_id) if dset_id else "last"
            dsets[dset_str] = data.Dataset(rundate, tech=tech, stage=stage, dataset_name=session, dataset_id=dset_id)

        # Call the writers
        plugins.call_one(package_name=__name__, plugin_name=writer, prefix=prefix, dset=dsets[dset_str])


def write_one(writer, dset, **writer_args):
    """Call one writer

    Args:
        writer (String):   Name of writer.
        dset (Dataset):    Model run data.
    """
    plugins.call_one(package_name=__name__, plugin_name=writer, dset=dset, **writer_args)
