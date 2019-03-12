"""Set up a plug-in architecture for data tables within a Dataset

To add a new table, simply create a new .py-file which defines a class that
inherits from the abstract class table.Table, and implements the missing
methods and properties.



"""

# System library imports
from datetime import datetime
import glob
import importlib
import json
import os

# Where imports
from where.data.dataset import Dataset
from where.data.table import Table
from where.lib import config
from where.lib import files
from where.lib import log


def list_datasets(rundate, tech, session, stage, **kwargs):
    """List datasets in a given dataset file

    Args:
        rundate:  Datetime, the model run date.
        tech:     String, the technique.
        stage:    String, the stage.
        kwargs:   Other arguments are passed to files.open.

    Returns:
        List of strings describing the datasets.
    """
    file_vars = dict(
        config.program_vars(rundate, tech, session=session, stage=stage, **kwargs), **config.date_vars(rundate)
    )

    try:
        with files.open("dataset_json", file_vars=file_vars) as fid:
            json_data = json.load(fid)
    except FileNotFoundError:
        return list()
        log.fatal(f"No data found for {tech.upper()} {stage} {rundate.strftime(config.FMT_date)}")

    return sorted(k for k in json_data.keys() if not k.startswith("_") and "/" in k)


def list_dataset_names_and_ids(rundate, tech, session, stage, **kwargs):
    """List names and ids of datasets in a given dataset file

    Args:
        rundate:  Datetime, the model run date.
        tech:     String, the technique.
        stage:    String, the stage.
        kwargs:   Other arguments are passed to files.open.

    Returns:
        2-tuple, list of strings with dataset names and list of ints with dataset ids.
    """
    datasets = list_datasets(rundate, tech, session, stage, **kwargs)
    if datasets:
        names, ids = zip(*[d.split("/") for d in datasets])
        return names, [int(i) for i in ids]
    else:
        return [], []


def list_dataset_names(rundate, tech, session, stage, **kwargs):
    """List names of datasets in a given dataset file

    Use list_names_and_ids and remove duplicates.

    Args:
        rundate:  Datetime, the model run date.
        tech:     String, the technique.
        stage:    String, the stage.
        kwargs:   Other arguments are passed to files.open.

    Returns:
        List of strings with dataset names.
    """
    names = list_dataset_names_and_ids(rundate, tech, session, stage, **kwargs)[0]
    return sorted(set(names))


def list_dataset_ids(rundate, tech, session, stage, name, **kwargs):
    """List ids of datasets with a given name in a given dataset file

    Use list_names_and_ids and filter on name.

    Args:
        rundate:  Datetime, the model run date.
        tech:     String, the technique.
        stage:    String, the stage.
        name:     String, the name of the dataset.
        kwargs:   Other arguments are passed to files.open.

    Returns:
        List of ints with dataset ids.
    """
    names, ids = list_dataset_names_and_ids(rundate, tech, session, stage, **kwargs)
    return [i for n, i in zip(names, ids) if n == name]


def list_datasets_from_filename(filename):
    """List all datasets from a given filename

    This is very hacky, and is not able to deal with user, sessionid, archived files etc.

    TODO: This does not work at the moment, because session is not handled
    """
    # Get info about dataset from filename
    tech, rundate_session, stage, _ = os.path.basename(filename).split(".")[0].split("-")

    session = rundate_session[8:]
    rundate = datetime.strptime(rundate_session[0:8], "%Y%m%d")
    names, ids = list_dataset_names_and_ids(rundate, tech, session, stage)

    dataset_list = list()
    for dataset_name, dataset_id in zip(names, ids):
        dataset_list.append(
            dict(
                rundate=rundate,
                tech=tech,
                session=session,
                stage=stage,
                dataset_name=dataset_name,
                dataset_id=dataset_id,
            )
        )
    return dataset_list


def parse_dataset_id(rundate, tech, stage, dataset_name, dataset_id, **kwargs):
    """Allow for some advanced handling of dataset_id

    In addition to using regular numbers as dataset_id, some text keywords can be used:

    + 'last': Use the last dataset_id written to file, default 0 if no file is previously written.
    + 'all':  Return a list of all dataset_ids in the file.
    """
    if isinstance(dataset_id, (float, int)):
        return dataset_id

    # Use the JSON-file to find information about the dataset ids
    file_vars = dict(
        config.program_vars(rundate, tech, session=dataset_name, stage=stage, **kwargs), **config.date_vars(rundate)
    )
    try:
        with files.open("dataset_json", file_vars=file_vars) as fid:
            json_data = json.load(fid)
    except FileNotFoundError:
        json_data = dict()

    if dataset_id == "last":
        # If _last_dataset_id is not given, use dataset_id=0 as default
        return json_data.get(dataset_name, dict()).get("_last_dataset_id", 0)

    if dataset_id == "all":
        return [int(k.split("/")[-1]) for k in json_data.keys() if k.startswith("{}/".format(dataset_name))]


def _register_tables():
    """Import data tables in the current directory and register them
    """
    directory = os.path.dirname(__file__)
    packagename = os.path.splitext(__name__)[0]
    for filename in glob.glob(os.path.join(directory, "*.py")):
        if filename.startswith("_") or filename.startswith("dataset3"):
            continue
        modulename = os.path.splitext(os.path.basename(filename))[0]
        module = importlib.import_module(packagename + "." + modulename)

        for clsname in dir(module):
            cls = getattr(module, clsname)
            try:
                is_table = issubclass(cls, Table) and cls.datatype is not None
            except TypeError:
                is_table = False
            if is_table:
                add_funcname = "add_" + cls.datatype
                read_funcname = "_read_" + cls.datatype
                log.debug(f"Registering data table {cls.datatype} from {filename} to Dataset.{add_funcname}")
                setattr(Dataset, add_funcname, Dataset._add(cls))
                setattr(Dataset, read_funcname, Dataset._read(cls))


# Register data tables to Dataset when module is loaded
_register_tables()
