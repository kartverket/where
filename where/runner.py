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

# Midgard imports
from midgard.dev.timer import Timer

# Where imports
import where
from where import pipelines
from where import setup
from where.lib import config
from where.lib import log
from where.lib import util
from where.lib.enums import LogLevel

_STATISTICS = {"Number of analyses": 0, "Successful analyses": 0, "Failed analyses": 0}


@Timer(f"Finish {util.get_program_name()} in")
def main():
    """Parse command line options and loop over the Where analysis

    Do simple parsing of command line arguments. Set up config-files and potentially start the analysis. See the help
    docstring at the top of the file for more information about the workflow.
    """
    util.check_help_and_version(doc_module=__name__)
    log.init(log_level=config.where.log.default_level.str, prefix="Runner")

    # Initialize
    pipeline = pipelines.get_from_options()
    config.read_pipeline(pipeline)
    if util.check_options("--doy"):
        from_date = util.parse_args("doy", doc_module=__name__)
        to_date = util.parse_args("doy", doc_module=__name__)
        sys.argv.remove("--doy")
    else:
        from_date = util.parse_args("date", doc_module=__name__)
        to_date = util.parse_args("date", doc_module=__name__)

    id_ = util.read_option_value("--id", default="")

    setup.set_profile(pipeline)

    # Initialize file variables
    file_vars = dict(
        pipeline=pipeline,
        id=id_,  
        **util.get_user_info(), 
        time=datetime.now().strftime("%Y%m%d%H%M%S"),
    )

    # Start logging
    log.file_init(
        file_path=config.files.path("log_runner", file_vars=file_vars),
        log_level=config.where.log.default_level.str,
        prefix="Runner",
        rotation=config.where.log.number_of_log_backups.int,
    )
    atexit.register(log_statistics)

    # Should where_runner crash if Where crashes?
    stop_on_error_opts = None
    if util.check_options("--stop-on-error"):
        stop_on_error_opts = True
    elif util.check_options("--continue-on-error"):
        stop_on_error_opts = False
    stop_on_error = config.where.get("stop_on_error", section="runner", value=stop_on_error_opts).bool
    error_logger = log.fatal if stop_on_error else log.error

    # Loop over dates
    rundate = from_date
    while rundate <= to_date:
        args = remove_runner_args(sys.argv[1:])
        where_args = set(pipelines.get_args(rundate, pipeline, input_args=args))

        for arg in sorted(where_args):
            cmd = f"{where.__executable__} {rundate:%Y %m %d} ".split() + arg.split()
            log.info(f"Running '{' '.join(cmd)}'")
            count("Number of analyses")
            try:
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as err:
                count("Failed analyses")
                error_msg = err.stderr.decode().strip().split("\n")[-1]
                error_logger(f"Command '{' '.join(cmd)}' failed: {error_msg}")
            else:
                count("Successful analyses")
            copy_log_from_where(rundate, pipeline, arg)

        rundate += timedelta(days=1)


def remove_runner_args(args):

    args = set(args)
    runner_args = set()
    for arg in args:
        option = arg.replace("-", "").split("=")[0]
        if config.where.exists(option, section="runner"):
            runner_args.add(arg)

    return list(args - runner_args)


def copy_log_from_where(rundate, pipeline, args):
    kwargs = dict()
    for a in args.split():
        if "=" in a:
            a = a.split("=", maxsplit=1)
            kwargs[a[0].lstrip("-")] = a[1]

    file_vars = dict(
        **config.program_vars(rundate, pipeline, use_options=False, **kwargs), **config.date_vars(rundate)
    )
    log_level = config.where.runner.log_level.str
    current_level = "none"
    try:
        with config.files.open("log", file_vars=file_vars) as fid:
            for line in fid:
                line_level, _, text = line.partition(" ")
                line_level = line_level.strip().lower()
                current_level = line_level if line_level else current_level
                text = text.strip()
                if getattr(LogLevel, current_level) >= getattr(LogLevel, log_level) and text:
                    # strip the 20 date characters from the text
                    log.log(text[20:], current_level)
    except FileNotFoundError as err:
        log.warn(f"'{err}'")


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
