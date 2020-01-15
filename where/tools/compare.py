#!/usr/bin/env python3
"""Compare Where datasets

Usage::

    {exe:tools} compare date --<pipeline> --items=<items> --specifier=<specifier> --stage=<stage> --writers=<writers> [options]

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
--items=             List with items to compare. Items can be 'id', 'session' or
                     'stage'. The given items are specified by 'specifier' 
                     option.
--specifier=         Specifier for given items, which can be 'id', 'session' or 
                     'stage'.
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
--id=                Special analysis identifier. Used in case that option 'ids'
                     is not defined.
-h, --help           Show this help message and exit.
===================  ===========================================================


Description:
------------

The script compares different Where datasets based on given writers. The 
different datasets are identified by the given items. The optione 'specifier'
defines, which kind of items are given.

Examples:
---------
Concatenate datasets from SISRE analysis by using 'id' items:
{exe:tools} compare 2019 1 1 2019 1 1 --sisre --stage=calculate --dset_id=1 --writers=sisre_comparison_report --items='grc_inav_e1_std_sat_concatenated,grc_inav_e1e5b_std_sat_concatenated,grc_fnav_e1e5a_std_sat_concatenated' --specifier=id

Concatenate datasets from SISRE analysis by using 'session' items:
{exe:tools} compare 2019 7 1 2019 7 1 --gnss --stage=estimate --writers=gnss_comparison_report --id=cnes_inav_e1_concatenated --items='nabf, hons, vegs, krss' --specifier=session

Concatenate datasets from SISRE analysis by using 'stage' items:
{exe:tools} compare 2019 1 1 2019 1 1 --gnss --dset_name=krss --dset_id=1 --writers=sisre_comparison_report --items='calculate, estimate' --specifier=stage

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
def main(date: "datedoy", tech: "pipeline", items: "option", specifier: "option"):
    log.init(log_level="info")
    dsets = dict()

    # Additional options
    stage = util.read_option_value("--stage")
    writer_names = util.read_option_value("--writers").replace(",", " ").split()
    items_ = [s.strip() for s in items.split(",")]

    # Get optional options
    dataset_id = util.read_option_value("--dset_id", default="last")
    dataset_id = "last" if dataset_id == "last" else int(dataset_id)
    dataset_name = util.read_option_value("--dset_name", default="")
    session = util.read_option_value("--session", default="")
    id_ = "-" + util.read_option_value("--id", default="") if util.read_option_value("--id", default="") else ""

    # Read datasets for given specifier
    if specifier == "id":
        for id_ in items_:
            dset = data.Dataset(
                rundate=date, tech=tech, stage=stage, dataset_name=dataset_name, dataset_id=dataset_id, id="-" + id_
            )
            if dset.num_obs == 0:
                log.warn(f"Dataset '{id_}' is empty.")
                continue
            dsets.update({id_: dset})

    elif specifier == "session":
        for session in items_:
            dset = data.Dataset(
                rundate=date, tech=tech, stage=stage, dataset_name=session, dataset_id=dataset_id, id=id_
            )
            if dset.num_obs == 0:
                log.warn(f"Dataset '{session}' is empty.")
                continue
            dsets.update({session: dset})

    elif specifier == "stage":
        for stage in items_:
            dset = data.Dataset(
                rundate=date, tech=tech, stage=stage, dataset_name=dataset_name, dataset_id=dataset_id, id=id_
            )
            if dset.num_obs == 0:
                log.warn(f"Dataset '{stage}' is empty.")
                continue
            dsets.update({stage: dset})
    else:
        log.fatal(f"Specifier {specifier} is not defined. It should be either 'id', 'session' or 'stage'.")

    if len(dsets) == 0:
        log.fatal(f"All given datasets are empty [{', '.join(dsets.keys())}].")
    elif len(dsets) == 1:
        log.warn(f"Nothing to compare. Only dataset '{list(dsets.keys())[0]}' is available.")

    # Loop over writers
    for writer in writer_names:
        write(writer, dset=dsets)


if __name__ == "__main__":
    main()
