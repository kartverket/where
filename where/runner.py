#!/usr/bin/env python3
"""Run the Where program to do analysis of space geodetic data over a time period

Usage::

    {exe:runner} from_date to_date --<pipeline> [options]

The program requires two dates, specifying the time period that will be
analysed. Typically, each date is given in the format `<year month day>` (for
example 2018 6 7). However, it is also possible to specify the date as `<year
day-of-year>` (for example 2015 216) by also adding the option `--doy`.  In
addition, one pipeline must be specified. See below for available pipelines.

The following pipelines are recognized:

===================  ===========================================================
Option               Description
===================  ===========================================================
{pipelines_doc:Run}

--doy                Specify from- and to-dates as Day-Of-Year
--stop-on-error      Stop runner if one analysis crashes.
--continue-on-error  Continue runner even if one analysis crashes.
--version            Show version information and exit.
-h, --help           Show this help message and exit.
===================  ===========================================================

In addition can all regular where options be used. See

    {exe:where} --help

for a list and explanation of those options.

Description:
------------

This program is used to run several Where analyses.


Examples:
---------

Run a VLBI analysis for all of 2015::

    {exe:runner} 2015 1 1 2015 12 31 -v


Current Maintainers:
--------------------

{maintainers}

Version: {version}

"""
# Standard library imports
import atexit
from datetime import datetime, timedelta
import subprocess
import sys


# Where imports
from midgard.dev.console import indent

import where
from where.lib import config
from where.lib import files
from where.lib import log
from where import pipelines
from where.lib.timer import timer
from where.lib import util


_STATISTICS = {"Number of analyses": 0, "Successful analyses": 0, "Failed analyses": 0}


@timer(f"Finish {util.get_program_name()} in")
def main():
    """Parse command line options and loop over the Where analysis

    Do simple parsing of command line arguments. Set up config-files and potentially start the analysis. See the help
    docstring at the top of the file for more information about the workflow.
    """
    # Initialize
    if util.check_options("--doy"):
        from_date = util.parse_args("doy", doc_module=__name__)
        to_date = util.parse_args("doy", doc_module=__name__)
    else:
        from_date = util.parse_args("date", doc_module=__name__)
        to_date = util.parse_args("date", doc_module=__name__)
    tech = pipelines.get_from_options()

    # Handle list of sessions
    session_list = set(util.read_option_value("--session", default="").replace(",", " ").split())
    sys.argv = [o for o in sys.argv if not o.startswith("--session=")]

    # Start logging
    log.init()
    file_vars = dict(timestamp=datetime.now().strftime(config.FMT_dt_file), **util.get_user_info())
    log.file_init(files.path("log_runner", file_vars=file_vars), config.where.runner.filelog_level.str)
    atexit.register(log_statistics)

    # Should where_runner crash if Where crashes?
    stop_on_error_opts = None
    if util.check_options("--stop-on-error"):
        stop_on_error_opts = True
    elif util.check_options("--continue-on-error"):
        stop_on_error_opts = False
    stop_on_error = config.where.get("stop_on_error", section="runner", value=stop_on_error_opts).bool
    error_logger = log.fatal if stop_on_error else log.error

    # The remaining options are passed onwhere to Where
    where_args = sys.argv[1:]

    # Loop over dates
    rundate = from_date
    while rundate <= to_date:
        available_sessions = set(pipelines.list_sessions(rundate, tech))
        sessions = available_sessions & session_list if session_list else available_sessions

        for session in sorted(sessions):
            cmd = f"{where.__executable__} {rundate:%Y %m %d} --session={session}".split() + where_args
            log.blank()
            log.blank(log_to_file=True)
            log.info(f"Running '{' '.join(cmd)}'")
            count("Number of analyses")
            try:
                subprocess.run(cmd, check=True, stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as err:
                count("Failed analyses")

                # Attempt to recover errors and traceback from stderr
                stderr = err.stderr.decode()
                stderr, *tb_exc = stderr.partition("\nTraceback")  # Regular Python traceback
                stderr, *tb_fatal = stderr.partition("\nFATAL:")  # Where log.fatal traceback
                stderr, *tb_error = stderr.partition("\nERROR:")  # Where log.error traceback
                traceback = indent("".join(tb_exc + tb_fatal + tb_error).strip(), 4)
                print(stderr, file=sys.stderr)
                error_logger(f"Command '{' '.join(cmd)}' failed with\n{traceback}")
            else:
                count("Successful analyses")
            copy_log_from_where(rundate, tech, session)

        rundate += timedelta(days=1)


def copy_log_from_where(rundate, tech, session):
    levels_to_log = config.where.runner.levels_to_log.list
    file_vars = dict(**config.program_vars(rundate, tech, session), **config.date_vars(rundate))
    log_path = files.path("log", file_vars=file_vars)
    if not log_path.exists():
        log.warn(f"Found no log at '{log_path}'")
        return

    log.blank()
    log.info(f"Information collected from '{log_path}':")
    with files.open("log", file_vars=file_vars) as fid:
        log_level = ""
        for line in fid:
            if not line.strip():  # Skip empty lines
                continue

            _log_level, _, message = line.rpartition(":")
            if _log_level in log.LOGLEVELS:
                log_level = _log_level.lower()
                message = message.strip()
                print_level = True
            else:  # Lines without a log level belong to the previous line
                message = ":".join((_log_level, message)).rstrip()
                print_level = False
            if log_level in levels_to_log:
                try:  # Temporary fix because .format sometimes complains
                    log.log(log_level.upper(), message, print_level=print_level)
                except KeyError:
                    pass


def count(statistic):
    _STATISTICS[statistic] += 1


def log_statistics():
    log.blank()
    log.info("Summary:")
    for statistic, value in _STATISTICS.items():
        log.info(f"{statistic}: {value}")


# Run main function only when running as script
if __name__ == "__main__":
    sys.exit(main())
