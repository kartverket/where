#!/usr/bin/env python3
"""Concatenate Where datasets

Usage::

    {exe:tools} concatenate <from_date> <to_date> --<pipeline> --stage=<stage> [options]

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
--stage=             Stage of analysis.
===================  ===========================================================

Furthermore, the following options are recognized:

===================  ===========================================================
Option               Description
===================  ===========================================================
--doy                Specify date as <year day-of-year>.
--dset_id=           Dataset identifier (Default: 'last').
--dset_name=         Dataset name (Default: '').
--id=                Analysis identifier (Default: '').
--only_for_rundate   Concatenate only data for given run date. Data are removed
                     from Datasets, which exceeds run date boundaries.
--session=           Session name (Default: '').
--writers=           List with writers.
-h, --help           Show this help message and exit.
===================  ===========================================================

Defined 'writers' uses configuration given in general Where configuration. If
needed each configuration option can be updated by adding command line options 
like:

    --<key>=<value>             # to update an entry in the 'pipeline' section
    --<section>:<key>=<value>   # to update an entry in a specific section

Example:
--------
Here are some concrete examples of how to run the program:

Concatenate datasets from GNSS analysis for a given period:
  {exe:tools} concatenate 2019 1 1 2019 1 31 --gnss --session=vegs --stage=estimate --id=grc_inav_e1 --writers=gnss_report

Concatenate datasets from SISRE analysis for a given period:
  {exe:tools} concatenate 2018 2 1 2018 2 3 --sisre --stage=calculate --id=mgex_fnav_e1e5a_5min --dset_id=2

Concatenate datasets from RINEX_NAV analysis for a given period:
  {exe:tools} concatenate 2018 2 1 2018 2 3 --rinex_nav --stage=edit --only_for_rundate --rinex_nav_report:only_for_rundate=False 
"""
# Standard library imports
from datetime import date, datetime, timedelta
from typing import Dict, List, Tuple
import sys

# External library imports
import numpy as np

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
def concatenate(from_date: "datedoy", to_date: "datedoy", tech: "pipeline", stage: "option"):
    log.init(log_level="info")

    # Get options
    dataset_id = util.read_option_value("--dset_id", default="last")
    dataset_id = "last" if dataset_id == "last" else int(dataset_id)
    dataset_name = util.read_option_value("--dset_name", default="")
    id = util.read_option_value("--id", default="")
    only_for_rundate = True if util.check_options("--only_for_rundate") else False
    session = util.read_option_value("--session", default="")
    writer_names = util.read_option_value("--writers", default="").replace(",", " ").split()

    # Update configuration of Where analysis
    config.where.update_from_options(_clean_sys_argv(tech))

    dset_vars = dict(
        tech=tech,
        stage=stage,
        session=session,
        dataset_name=session,
        dataset_id=dataset_id,
        session_name=id + "_concatenated",
    )
    dset = _concatenate_datasets(from_date, to_date, dset_vars, only_for_rundate)
    if dset.num_obs == 0:
        log.fatal(f"No data to read period from {from_date} to {to_date}.")
    dset.write()

    # Loop over writers
    for writer in writer_names:
        write(writer, dset=dset)


#
# AUXILIARY FUNCTIONS
#
def _clean_sys_argv(tech: str) -> List[str]:
    """Values in sys.argv that are not valid option values in Where
    """
    reserved_opts = {tech, "dset_id", "dset_name", "id", "only_for_rundate", "session", "stage", "writers"}
    return [o for o in sys.argv[1:] if o.startswith("--") and o[2:].split("=")[0] not in reserved_opts]


def _concatenate_datasets(
    from_date: date, to_date: date, dset_vars: Dict[str, str], only_for_rundate: bool
) -> np.ndarray:
    """Concatenate datasets

    Args:
        from_date:         Start date for reading Dataset.
        to_date:           End date for reading Dataset.
        dset_vars:         Common Dataset variables.
        only_for_rundate:  Concatenate only data for given rundate.
    """
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

        if only_for_rundate:
            _keep_data_only_for_rundate(dset)

            if dset.num_obs == 0:
                log.info(f"No data to read for {current_date}")
                continue

        log.info(f"Reading data for {current_date}")
        if not dset_merged:
            dset_merged.copy_from(dset)
        else:
            dset_merged.extend(dset)

    return dset_merged


def _get_day_limits(rundate: date) -> Tuple[datetime, datetime]:
    """Get start and end time for given run date

    Args:
        rundate: Run date.

    Returns:
        Start and end date. 
    """
    day_start = datetime(rundate.year, rundate.month, rundate.day)
    day_end = day_start + timedelta(seconds=86399.9999)

    return day_start, day_end


def _keep_data_only_for_rundate(dset: "Dataset") -> None:
    """Keep only data of Dataset from the current run date

    Args:
        dset: A dataset containing the data.
    """
    day_start, day_end = _get_day_limits(dset.analysis["rundate"])

    # Remove data not belonging to current rundate
    idx = np.logical_and(dset.time.datetime >= day_start, dset.time.datetime <= day_end)
    dset.subset(idx)
