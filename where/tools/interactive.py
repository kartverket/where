"""Do an interactive Where analysis

Description:
------------

The :func:`interactive` function reads data from previous model runs and makes
them available in an interactive IPython session.

"""
# External library imports
import IPython

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import data
from where.lib import config
from where.lib import files

# Other imports useful for the interactive session
import numpy as np  # noqa
import pandas as pd  # noqa
import matplotlib.pyplot as plt  # noqa
from where import apriori  # noqa
from midgard.math.constant import constant  # noqa
from where.lib.time import Time, TimeDelta  # noqa
from where.lib.unit import Unit  # noqa


@plugins.register
def interactive(rundate: "date", pipeline: "pipeline", session: "option" = ""):  # typing: ignore
    """Read model run data and start an interactive session

    Read all datasets for the given rundate and pipelines, and start an interactive IPython session.

    Args:
        rundate: The model run date.
        pipeline:    String with the name of pipelinenique.
    """
    # Read data for all pipelines
    vars_dict = dict()
    list_of_vars = list()
    config.init(rundate, pipeline, session=session)
    vars_dict[pipeline] = config.files.vars.copy()
    del vars_dict[pipeline]["rundate"]

    # Register filekey suffix
    filekey_suffix = config.tech.get("filekey_suffix", default="").list
    if filekey_suffix:
        config.files.profiles = filekey_suffix

    # Read data for all available sessions and stages, add them to the global namespace
    stages = files.glob_variable("dataset_hdf5", "stage", ".+")
    for stage in sorted(stages):
        names, dset_ids = data.list_dataset_names_and_ids(rundate, stage=stage, **vars_dict[pipeline])
        for name, dset_id in zip(names, dset_ids):
            var_name = "_".join([pipeline, stage, name, str(dset_id)])
            short_var_name = pipeline[0] + str(len([v for v in list_of_vars if v.lstrip().startswith(pipeline[0])]))
            globals()[var_name] = data.Dataset(
                rundate, stage=stage, dataset_name=name, dataset_id=dset_id, **vars_dict[pipeline]
            )
            globals()[short_var_name] = globals()[var_name]
            list_of_vars.append("{:>6s}, {}".format(short_var_name, var_name))

    # Start an interactive IPython session
    IPython.embed(header="Available datasets:\n" + "\n".join(list_of_vars))
