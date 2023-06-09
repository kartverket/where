#!/usr/bin/env python3
"""Run the Where program to do analysis of space geodetic data

Usage:

    {exe} <date> <pipeline> [options]

The program requires a date. Typically, the date is given in the format
`<year month day>` (for example 2015 8 4). However, it is also possible to
specify the date as `<year day-of-year>` (for example 2015 216) by also adding
the option `--doy`.

In addition, one pipeline must be specified. See below for available pipelines.

===================  ===========================================================
Pipeline             Description
===================  ===========================================================
{pipelines_doc:Run}
===================  ===========================================================

Furthermore, the following options are recognized:

===================  ===========================================================
Option               Description
===================  ===========================================================
-D, --delete         Delete existing analysis results.
    --doy            Specify date as <year day-of-year>.
-E, --edit           Edit the configuration of an analysis.
-F, --force          Run all stages even if no dependencies have changed.
-I, --interactive    Start an interactive session with analysis data available.
-N, --new            Start an analysis with a new config.
-S, --showconfig     Show the configuration of an analysis.
-T, --showtb         Show traceback if the program crashes.
--id=analysisid      Add a special analysis id (to run several versions of the
                     same analysis simultaneously).
--profile=name       Use config settings specfied for a given profile, for
                     instance --profile=vascc for a VLBI analysis.
--session_code=session_code    Run analysis for the given session.
--station=station    Run analysis for the given station. This is only valid for
                     GNSS, GNSS_OBS and GNSS_VEL pipeline.
--user=username      Run as username. Does not need to be an existing username
                     on the system.
--debug, ...         Show additional debug information. Other flags such as
                     --all, --debug, --time, --dev, --info, --out, --warn,
                     --check, --error, --fatal, --none are also allowed, and
                     shows differing amounts of information as the program runs.
--version            Show version information and exit.
-h, --help           Show this help message and exit.
===================  ===========================================================

Finally, configuration settings of an analysis can be changed using command line
options. Run an analysis with the `-S` option for details.


Description:
------------

This program is used to run a Where analysis.


Examples:
---------

Here are some concrete examples of how to run Where:

Run a VLBI analysis for August 4 2015:

    {exe} 2015 8 4 -v --session_code=R1699

Run an SLR analysis for September 1 2015 using day-of-year:

    {exe} 2015 242 --slr --doy

Change the sampling_rate option of the GNSS analysis:

    {exe} 2016 3 1 -g --station=stas --sampling_rate=30

Look at the configuration of VLBI analysis set up for November 2 2009:

    {exe} 2009 11 2 --vlbi -S --session_code=R1403


Current Maintainers:
--------------------

{maintainers}

Version: {version}

"""
# Standard library imports
import sys

# Midgard imports
from midgard.dev.timer import Timer

# Where imports
from where import pipelines
from where import setup
from where.lib import config
from where.lib import log
from where.lib import util


@Timer(f"Finish {util.get_program_name()} in")
@util.no_traceback
def main():
    """Parse command line options and run the Where analysis

    Do simple parsing of command line arguments. Set up config-files and start the analysis. See the help docstring at
    the top of the file for more information about the workflow.
    """
    util.check_help_and_version(doc_module=__name__)

    # Start logging
    log.init(config.where.log.default_level.str)
    log.debug(f"Use {util.get_python_version()} on process {util.get_pid_and_server()}")

    # Read command line options
    pipeline = pipelines.get_from_options()
    config.read_pipeline(pipeline)
    if util.check_options("--doy"):
        rundate = util.parse_args("doy", doc_module=__name__)
    else:
        rundate = util.parse_args("date", doc_module=__name__)

    args, kwargs = util.options2args(sys.argv[1:])

    # Start an interactive session
    if util.check_options("-I", "--interactive"):
        from where.tools import interactive  # Local import because interactive imports many external packages

        interactive.interactive(rundate, pipeline, **kwargs)
        return

    # Set up the configuration for a new analysis or update an existing one
    unused_options = setup.setup_config(rundate, pipeline, *args, **kwargs)

    pipeline_args, pipeline_kwargs = util.options2args(unused_options)

    # Run the analysis
    setup.add_timestamp(rundate, pipeline, "last run", **kwargs)
    with Timer(f"Finish pipeline {pipeline.upper()} in"):
        pipelines.run(rundate, pipeline, *pipeline_args, **pipeline_kwargs)


# Run main function only when running as script
if __name__ == "__main__":
    sys.exit(main())
