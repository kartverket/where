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

# Midgard imports
from midgard.dev import plugins
from midgard import writers as mg_writers
from midgard.writers import names, write as write_one  # noqa

# Where imports
from where.lib import config
from where import data

# Add Where writers to Midgard writers
plugins.add_alias(mg_writers.__name__, __name__)


def write(default_dset):
    """Call all writers specified in the configuration

    The list of writers to use is taken from the config file of the given technique. Each writer is passed a
    :class:`~where.data.dataset.Dataset` with data for the modelrun and should write the relevant parts of the data to
    file.

    By default the last dataset for the default_stage is sent to the writer, but that is possible to override with the
    following notation:

        output = writer_1                 # Use default dataset
        output = writer_1:calculate       # Use last dataset of "calculate" stage
        output = writer_1:calculate/2     # Use dataset 2 of "calculate" stage

    Args:
        default_dset (Dataset):    Dataset used by default.
    """
    dsets = {f"{default_dset.vars['stage']}/{default_dset.vars['dataset_id']}": default_dset}
    prefix = config.analysis.get("analysis", default="").str
    output_list = config.tech.output.list
    writer_and_dset = [o.partition(":")[::2] for o in output_list]

    rundate = config.analysis.rundate.date
    tech = config.analysis.tech.str
    session = config.analysis.session.str
    for writer, dset_str in writer_and_dset:

        # Read the datasets
        if dset_str not in dsets:
            stage, _, dset_id = dset_str.partition("/")
            stage, _, dset_name = stage.partition(":")
            stage = stage if stage else default_dset.vars["stage"]
            dset_name = dset_name if dset_name else session
            dset_id = int(dset_id) if dset_id else "last"
            dsets[dset_str] = data.Dataset(
                rundate, tech=tech, stage=stage, dataset_name=dset_name, dataset_id=dset_id, session=session
            )

        # Call the writers
        plugins.call(package_name=mg_writers.__name__, plugin_name=writer, prefix=prefix, dset=dsets[dset_str])
