#!/usr/bin/env python3
"""Do an interactive Where analysis

Description:
------------

The :func:`interactive` function reads data from previous model runs and makes
them available in an interactive IPython session.

"""

# External library imports
import IPython

# Where imports
from where import data
from where.lib import config
from where.lib import files

# Other imports useful for the interactive session
import numpy as np  # noqa
import pandas as pd  # noqa
import matplotlib.pyplot as plt  # noqa
from where import apriori  # noqa
from where.lib import constant  # noqa
from where.lib.time import Time, TimeDelta  # noqa
from where.lib.unit import unit  # noqa


def interactive(rundate, tech, session=""):
    """Read model run data and start an interactive session

    Read all datasets for the given rundate and techniques, and start an interactive IPython session.

    Args:
        rundate: The model run date.
        tech:    String with the name of technique.
    """
    # Read data for all techniques
    vars_dict = dict()
    list_of_vars = list()
    config.init(rundate, tech, session=session)
    vars_dict[tech] = config.files.vars.copy()
    del vars_dict[tech]["rundate"]

    # Register filekey suffix
    filekey_suffix = config.tech.get("filekey_suffix", default="").list
    if filekey_suffix:
        files.use_filelist_profiles(*filekey_suffix)

    # Read data for all available sessions and stages, add them to the global namespace
    stages = files.glob_variable("dataset_hdf5", "stage", ".+")
    for stage in sorted(stages):
        names, dset_ids = data.list_dataset_names_and_ids(rundate, stage=stage, **vars_dict[tech])
        for name, dset_id in zip(names, dset_ids):
            var_name = "_".join([tech, stage, name, str(dset_id)])
            short_var_name = tech[0] + str(len([v for v in list_of_vars if v.lstrip().startswith(tech[0])]))
            globals()[var_name] = data.Dataset(
                rundate, stage=stage, dataset_name=name, dataset_id=dset_id, **vars_dict[tech]
            )
            globals()[short_var_name] = globals()[var_name]
            list_of_vars.append("{:>6s}, {}".format(short_var_name, var_name))

    # Start an interactive IPython session
    IPython.embed(header="Available datasets:\n" + "\n".join(list_of_vars))
