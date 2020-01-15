"""Framework for writing output in different formats

Description:
------------

Each output format / output destination should be defined in a separate .py-file. The function inside the .py-file that
should be called need to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midgard.dev import plugins

    @plugins.register
    def write_as_fancy_format(dset):
        ...

The decorated function will be called with a single parameter, ``dset`` which contains a
:class:`~where.data.dataset.Dataset` with the data that should be output.

"""

# Standard library imports
import os
import pathlib
import shutil

# Midgard imports
from midgard.dev import plugins
from midgard import writers as mg_writers
from midgard.writers import names, write as write_one  # noqa

# Where imports
from where.lib import config
from where.lib import log
from where.data import dataset3 as dataset

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
    prefix = config.analysis.get("analysis", default="").str
    output_list = config.tech.output.list
    writer_and_dset = [o.partition(":")[::2] for o in output_list]

    dset_vars = config.analysis.config.as_dict()
    dset_vars["rundate"] = config.analysis.rundate.date

    for writer, dset_str in writer_and_dset:
        # Read the datasets
        if dset_str:
            stage, _, label = dset_str.partition("/")
            stage = stage if stage else default_dset.vars["stage"]
            label = label if label else "last"
            dset = dataset.Dataset.read(stage=stage, label=label, **dset_vars)
            plugins.call(package_name=mg_writers.__name__, plugin_name=writer, prefix=prefix, dset=dset)
        else:
            plugins.call(package_name=mg_writers.__name__, plugin_name=writer, prefix=prefix, dset=default_dset)

    publish_files()


def publish_files(publish=None):
    """Publish files to specified directories

    The publish string should list file_keys specified in files.conf. Each file_key needs to have a field named publish
    specifying a directory the file should be copied to.

    Args:
        publish (String):   List of file_keys that will be published.
    """
    if not config.where.files.publish.bool:
        return

    publish_list = config.tech.get("files_to_publish", value=publish).list
    for file_key in publish_list:
        try:
            source = config.files.path(file_key)
        except KeyError:
            log.error(f"File key '{file_key}' in publish configuration is unknown. Ignored")
            continue
        if not source.exists():
            try:
                log.error(f"File '{source}' (file key='{file_key}') does not exist, and can not be published")
            except KeyError:
                log.error(f"File key='{file_key}' has incomplete filename information and can not be published")
            continue

        try:
            destinations = config.files[file_key].publish.replaced.as_list(convert=pathlib.Path)
        except AttributeError:
            log.error(f"File key '{file_key}' does not specify 'publish' directory in file configuration. Ignored")
            continue

        # Copy file to destinations
        for destination in destinations:
            destination_path = destination / source.name
            try:
                if destination_path.exists():
                    os.remove(destination_path)
                destination.mkdir(parents=True, exist_ok=True)
                shutil.copy(source, destination)
            except (OSError, PermissionError, FileNotFoundError) as err:
                log.error(f"Unable to publish {file_key}-file {source} to {destination} because of {err}")
            else:
                log.info(f"Published {file_key}-file {source} to {destination}")
