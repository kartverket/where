#!/usr/bin/env python3
"""Differentiate data between two given Datasets

Usage::

    {exe:tools} difference date --<pipeline> --ids=<identifiers> --stage=<stage> --writers=<writers> [options]

The program requires dates. Typically, the date is given in the format
`<year month day>` (for example 2015 8 4). However, it is also possible to
specify the date as `<year day-of-year>` (for example 2015 216) by also adding
the option `--doy`.

The following commands are required:

===================  ===========================================================
Command              Description
===================  ===========================================================
date                 The model run date in the format ``<year month day>``.
{pipelines_doc:Differenciate results from}
--ids=               List with special analysis identifiers.
--difference_by=     Name of fields to be used for differentiation of the two 
                     Datasets (e.g. 'time, satellite').
--stage=             Stage of analysis.
===================  ===========================================================

Furthermore, the following options are recognized:

===================  ===========================================================
Option               Description
===================  ===========================================================
--doy                Specify date as <year day-of-year>.
--label=             Dataset label (Default: 'last').
--station=           Station name (Default: '').
--writers=           List with writers.
-h, --help           Show this help message and exit.
===================  ===========================================================


Description:
------------

The script differentiate two given Where datasets. The two datasets are identified 
by the given identifiers via option --ids.

Examples:
---------
{exe:tools} difference 2018 2 1 --sisre --stage=calculate --dset_id=2 --ids='inav, fnav' --difference_by='time, satellite'

"""
# Midgard imports
from midgard.dev import plugins
from midgard.writers import write

# Where imports
from where import data
from where.lib import log
from where.lib import util


@plugins.register
def main(date: "datedoy", pipeline: "pipeline", ids: "option"):
    log.init(log_level="info")

    # Additional required options
    identifiers = [id_.strip() for id_ in ids.split(",")]
    difference_by = util.read_option_value("--difference_by").replace(",", " ").split()
    stage = util.read_option_value("--stage")

    # Get optional options
    dataset_id = util.read_option_value("--dset_id", default="last")
    dataset_id = "last" if dataset_id == "last" else int(dataset_id)
    dataset_name = util.read_option_value("--dset_name", default="")
    writer_names = util.read_option_value("--writers", default="").replace(",", " ").split()
    station = util.read_option_value("--station", default="")

    # Get datasets
    dset = data.Dataset(
            rundate=date, 
            pipeline=pipeline, 
            stage=stage, 
            station=station, 
            label=label, 
            id="-" + identifiers[0],
    )

    dset_other = data.Dataset(
            rundate=date, 
            pipeline=pipeline, 
            stage=stage, 
            station=station, 
            label=label, 
            id="-" + identifiers[1],
    )

    if dset.num_obs == 0:
        log.warn(f"Nothing to differentiate. Dataset '{identifiers[0]}' is empty.")
        return 1

    if dset_other.num_obs == 0:
        log.warn(f"Nothing to differentiate. Dataset '{identifiers[1]}' is empty.")
        return 1

    # Differentiate dataset
    dset_diff = dset.difference(dset_other, index_by=','.join(difference_by))
    dset_diff.write_as(stage="difference")

    # Loop over writers
    for writer in writer_names:
        write(writer, dset_diff)


if __name__ == "__main__":
    main()
