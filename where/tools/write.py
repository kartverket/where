#!/usr/bin/env python3
"""Write Where dataset fields by selecting the used "writers"

Usage::

    {exe:tools} write <date> <pipeline> --stage=<stage> --writers=<writers> [options]

The program requires a date. Typically, the date is given in the format
`<year month day>` (for example 2015 8 4). However, it is also possible to
specify the date as `<year day-of-year>` (for example 2015 216) by also adding
the option `--doy`.

The following commands are required:

===================  ===========================================================
Command              Description
===================  ===========================================================
date                 The starting date in the format ``<year month day>``.
{pipelines_doc:Write results from}
--stage=             Stage of analysis.
--writers=           List with writers.
===================  ===========================================================

Furthermore, the following options are recognized:

===================  ===========================================================
Option               Description
===================  ===========================================================
--doy                Specify date as <year day-of-year>.
--label=             Dataset label (Default: 'last').
--id=                Analysis identifier (Default: '').
--station=           Station name (Default: '').
-h, --help           Show this help message and exit.
===================  ===========================================================

Description:
------------

Read Where dataset and write Where datasets fields via a defined writer.

Example:
--------
  {exe:tools} write 2021 7 1 --sisre --stage=calculate --label=raw --id='_cnes_inav_e1' --writers=sisre_plot
  {exe:tools} write 2019 2 1 --gnss --station=stas --stage=write --writers=gnss_report
  {exe:tools} write 2019 2 1 --rinex_nav --stage=read --writers=rinex_nav_report
"""
# Standard library imports
from typing import List
import sys

# Midgard imports
from midgard.dev import plugins
from midgard.writers import write as write_

# Where imports
from where.data import dataset3 as dataset
from where.lib import config
from where.lib import log
from where.lib import util


@plugins.register
def write(rundate: "datedoy", pipeline: "pipeline", stage: "option", writers: "option"):
    log.init(log_level="info")

    # Get options
    label = util.read_option_value("--label", default="None")
    # TODO: label = "last" if label == "last" else label
    id_ = util.read_option_value("--id", default="")
    station = util.read_option_value("--station", default="")
    writers = writers.replace(",", " ").split()

    # Update configuration of Where analysis
    config.where.update_from_options(_clean_sys_argv(pipeline))

    dset_vars = dict(pipeline=pipeline, stage=stage, station=station, label=label, id=id_)

    try:
        dset = dataset.Dataset.read(**dict(dset_vars, rundate=rundate))
        path = config.files.path("dataset", file_vars={**dset.vars, **dset.analysis})
        log.info(f"Read Where dataset files {path}.")
    except OSError as err:
        log.fatal(f"Unable to read data for {rundate}: {err}")

    # Loop over writers
    for writer in writers:
        log.info(f"Apply writer '{writer}'.")
        write_(writer, dset=dset)


#
# AUXILIARY FUNCTIONS
#
def _clean_sys_argv(pipeline: str) -> List[str]:
    """Values in sys.argv that are not valid option values in Where
    """
    reserved_opts = {pipeline, "label", "id", "stage", "station", "writers"}
    return [o for o in sys.argv[1:] if o.startswith("--") and o[2:].split("=")[0] not in reserved_opts]
