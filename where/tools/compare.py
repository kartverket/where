#!/usr/bin/env python3
"""Compare Where datasets

Usage::

    {exe:tools} compare date --<pipeline> --ids=<identifiers> --stage=<stage> --writers=<writers> [options]

The program requires dates. Typically, the date is given in the format
`<year month day>` (for example 2015 8 4). However, it is also possible to
specify the date as `<year day-of-year>` (for example 2015 216) by also adding
the option `--doy`.

The following commands are required:

===================  ===========================================================
Command              Description
===================  ===========================================================
date                 The model run date in the format ``<year month day>``.
{pipelines_doc:Plot results from}
--ids=               List with special analysis identifiers.
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


Description:
------------

The script compares different Where datasets based on given writers. The 
different datasets are identified by the given identifiers via option --ids.

Examples:
---------
{exe:tools} compare 2018 2 1 --sisre --stage=calculate --writers=sisre_comparison_report --dset_id=2 --ids='mgex_inav_e1e5b_5min_concatenated, mgex_fnav_e1e5a_5min_concatenated'

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
def main(date: "datedoy", tech: "pipeline", ids: "option"):
    log.init(log_level="info")

    # Additional required options
    stage = util.read_option_value("--stage")
    writer_names = util.read_option_value("--writers").replace(",", " ").split()
    identifiers = [id_.strip() for id_ in ids.split(",")]

    # Get optional options
    dataset_id = util.read_option_value("--dset_id", default="last")
    dataset_id = "last" if dataset_id == "last" else int(dataset_id)
    dataset_name = util.read_option_value("--dset_name", default="")
    session = util.read_option_value("--session", default="")

    # Loop over different dataset identifiers
    dsets = dict()
    for id_ in identifiers:
        dset = data.Dataset(
            rundate=date, tech=tech, stage=stage, dataset_name=dataset_name, dataset_id=dataset_id, id="-" + id_
        )
        if dset.num_obs == 0:
            log.warn(f"Dataset '{id_}' is empty.")
            continue
        dsets.update({id_: dset})

    if len(dsets) == 0:
        log.fatal(f"All given datasets are empty [{', '.join(dsets.keys())}].")
    elif len(dsets) == 1:
        log.warn(f"Nothing to compare. Only dataset '{list(dsets.keys())[0]}' is available.")

    # Loop over writers
    for writer in writer_names:
        write(writer, dset=dsets)


if __name__ == "__main__":
    main()
