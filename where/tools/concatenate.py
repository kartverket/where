#!/usr/bin/env python3
"""Concatenate Where datasets

Usage::

    {exe:tools} concatenate <from_date> <to_date> --<pipeline> --id=<identifier> --stage=<stage> --writers=<writers> [options]

The program requires dates. Typically, the date is given in the format
`<year month day>` (for example 2015 8 4). However, it is also possible to
specify the date as `<year day-of-year>` (for example 2015 216) by also adding
the option `--doy`.


The following commands are required:

===================  ===========================================================
Command              Description
===================  ===========================================================
from_date            The starting date in the format ``<year month day>``.
to_date              The ending date in the format ``<year month day>``.
{pipelines_doc:Plot results from}
--id=                Analysis identifier.
--stage=             Stage of analysis.
--writers=           List with writers.
===================  ===========================================================

Furthermore, the following options are recognized:

===================  ===========================================================
Option               Description
===================  ===========================================================
--doy                Specify date as <year day-of-year>.
--dset_id=           Dataset identifier (Default: 'last').
--dset_name=         Dataset name (Default: '').
--session=           Session name (Default: '').
-h, --help           Show this help message and exit.
===================  ===========================================================

Example:
--------

  {exe:tools} concatenate 2018 2 1 2018 2 3 --sisre --stage=calculate --id=mgex_fnav_e1e5a_5min --dset_id=2
"""
# Standard library imports
from datetime import timedelta

# Midgard imports
from midgard.dev import plugins
from midgard.writers import write

# Where imports
from where import pipelines
from where.lib import config
from where import data
from where.lib import log
from where.lib import util


@plugins.register
def concatenate(
    from_date: "datedoy", to_date: "datedoy", tech: "pipeline", stage: "option", writers: "option", id: "option"
):
    log.init(log_level="info")

    # Get options
    writer_names = writers.replace(",", " ").split()
    dataset_id = util.read_option_value("--dset_id", default="last")
    dataset_id = "last" if dataset_id == "last" else int(dataset_id)
    dataset_name = util.read_option_value("--dset_name", default="")
    session = util.read_option_value("--session", default="")

    dset_vars = dict(
        tech=tech,
        stage=stage,
        session=session,
        dataset_name=session,
        dataset_id=dataset_id,
        session_name=id + "_concatenated",
    )
    dset = concatenate_datasets(from_date, to_date, dset_vars)
    if dset.num_obs == 0:
        log.fatal(f"No data to read period from {from_date} to {to_date}.")
    dset.write()

    # Loop over writers
    for writer in writer_names:
        write(writer, dset=dset)


def concatenate_datasets(from_date, to_date, dset_vars):
    merged_vars = config.program_vars(rundate=from_date, tech_name=dset_vars["tech"], **dset_vars)
    merged_vars["id"] += "_concatenated"
    dset_merged = data.Dataset(**dict(merged_vars, rundate=from_date, empty=True))

    date_to_read = from_date
    while date_to_read <= to_date:
        dset = data.Dataset(rundate=date_to_read, **dset_vars)
        current_date = date_to_read
        date_to_read += timedelta(days=1)
        if dset.num_obs == 0:
            log.info(f"No data to read for {current_date}")
            continue
        log.info(f"Reading data for {current_date}")
        if not dset_merged:
            dset_merged.copy_from(dset)
        else:
            dset_merged.extend(dset)

    return dset_merged
