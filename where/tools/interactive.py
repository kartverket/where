"""Do an interactive Where analysis

Description:
------------

The :func:`interactive` function reads data from previous model runs and makes
them available in an interactive IPython session.

"""
# External library imports
import IPython

# Standard library imports
import sys

# Midgard imports
from midgard.dev import plugins
from midgard.math.constant import constant  # noqa

# Where imports
from where.data import dataset3 as dataset
from where.lib import config

# Other imports useful for the interactive session
import numpy as np  # noqa
import pandas as pd  # noqa
import matplotlib.pyplot as plt  # noqa
from where import apriori  # noqa
from where.data.time import Time, TimeDelta  # noqa
from where.lib.unit import Unit  # noqa


@plugins.register
def interactive(rundate: "date", pipeline: "pipeline", **kwargs):  # typing: ignore
    """Read model run data and start an interactive session

    Read all datasets for the given rundate and pipelines, and start an interactive IPython session.

    Args:
        rundate:     The model run date.
        pipeline:    String with the name of pipeline
    """
    # Read data for all pipelines
    vars_dict = dict()
    datasets = list()
    config.init(rundate, pipeline, **kwargs)
    vars_dict[pipeline] = config.files.vars.copy()
    del vars_dict[pipeline]["rundate"]
    del vars_dict[pipeline]["pipeline"]

    # Register filekey suffix
    filekey_suffix = config.tech.get("filekey_suffix", default="").list
    if filekey_suffix:
        config.files.profiles = filekey_suffix

    # Read data for all available sessions and stages, add them to the global namespace
    stages = config.files.glob_variable("dataset", "stage", ".+")
    for stage in sorted(stages):
        file_vars = vars_dict[pipeline].copy()
        file_vars.update(stage=stage)
        labels = config.files.glob_variable("dataset", "label", ".+", file_vars)
        for label in labels:
            file_vars["label"] = label
            var_name = config.files.path("dataset", file_vars=file_vars)
            short_var_name = pipeline[0] + str(len([v for v in datasets if v.lstrip().startswith(pipeline[0])]))
            globals()[var_name] = dataset.Dataset.read(
                rundate=rundate, pipeline=pipeline, stage=stage, label=label, **vars_dict[pipeline]
            )
            globals()[short_var_name] = globals()[var_name]
            datasets.append("{:>6s}, {}".format(short_var_name, var_name))

    # Start an interactive IPython session
    IPython.embed(header="Available datasets:\n" + "\n".join(datasets))
