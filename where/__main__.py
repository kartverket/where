#!/usr/bin/env python3
"""Run the Where program to do analysis of space geodetic data

Usage::

    {exe} date pipeline [--session=session] [options]

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
-A, --archive        Move existing analysis results to archive.
-D, --delete         Delete existing analysis results.
    --doy            Specify date as <year day-of-year>.
-E, --edit           Edit the configuration of an analysis.
-F, --force          Run all stages even if no dependencies have changed.
-I, --interactive    Start an interactive session with analysis data available.
-N, --new            Start a new analysis (combine with -A or -D).
-S, --showconfig     Show the configuration of an analysis.
-T, --showtb         Show traceback if the program crashes.
--id=analysisid      Add a special analysis id (to run several versions of the
                     same analysis simultaneously).
--profile=name       Use config settings specfied for a given profile, for
                     instance --profile=vascc for a VLBI analysis.
--session=session    Run analysis for the given session.
--user=username      Run as username. Does not need to be an existing username
                     on the system.
--debug, ...         Show additional debug information. Other flags such as
                     --all, --debug, --time, --dev, --info, --warn, --check,
                     --error, --fatal, --none are also allowed, and will show
                     differing amounts of information as the program runs.
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

Here are some concrete examples of how to run the program:

Run a VLBI analysis for August 4 2015:

    {exe} 2015 8 4 -v --session=XA

Run an SLR analysis for September 1 2015 using day-of-year:

    {exe} 2015 242 --slr --doy

Change the spam option of the GNSS analysis:

    {exe} 2016 3 1 -g --spam=ham

Look at the configuration of all VLBI analyses set up for November 2 2009::

    {exe} 2009 11 2 --vlbi -S

Run a complete SLR analysis for January 28 2017::

    {exe} 2017 1 28 -s -F

Run the VLBI analysis for August 4 2015 only for session XA::

    {exe} 2015 8 4 --vlbi --only_session=XA


Current Maintainers:
--------------------

{maintainers}

Version: {version}

"""
# Standard library imports
import sys

# Where imports
from where import pipelines
from where import setup
from where.lib import log
from where.lib.timer import timer
from where.lib import util


@timer(f"Finish {util.get_program_name()} in")
@util.no_traceback
def main():
    """Parse command line options and run the Where analysis

    Do simple parsing of command line arguments. Set up config-files and start the analysis. See the help docstring at
    the top of the file for more information about the workflow.
    """
    # Start logging
    log.init()

    # Read command line options
    if util.check_options("--doy"):
        rundate = util.parse_args("doy", doc_module=__name__)
    else:
        rundate = util.parse_args("date", doc_module=__name__)
    pipeline = pipelines.get_from_options()
    session = util.read_option_value("--session", default="")

    # Pretend to empty mailbox
    pretend_to_empty_mailbox()

    # Start an interactive session
    if util.check_options("-I", "--interactive"):
        from where.tools import interactive

        interactive.interactive(rundate, pipeline, session)
        return

    # Set up the configuration for a new analysis or update an existing one
    setup.setup_config(rundate, pipeline, session)

    # Run the analysis
    setup.add_timestamp(rundate, pipeline, session, "last run")
    with timer(f"Finish pipeline {pipeline.upper()} in"):
        pipelines.run(rundate, pipeline, session)


def pretend_to_empty_mailbox():
    """Pretend to try to empty a mailbox

    This function simply prints an error message to stderr. No actual attempt to empty a mailbox is done. This is
    included purely for nostalgic (and possibly some misunderstood backwards compatibility) reasons :P
    """
    if util.check_options("-M", "--mailbox"):
        print("/var/spool/mail/pha: Permission denied.", file=sys.stderr)


# Run main function only when running as script
if __name__ == "__main__":
    sys.exit(main())
