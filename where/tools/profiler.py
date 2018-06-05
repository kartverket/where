#!/usr/bin/env python3
"""Run the Where program to do analysis of space geodetic data

Usage::

    {exe} date pipeline [options]

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
--id=sessionid       Add a special session id (to run several analyses for the
                     same day simultaneously).
--only_session=s     Only run session <s>.
--profile=name       Use config settings specfied for a given profile, for
                     instance --profile=vascc for a VLBI analysis.
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

    {exe} 2015 8 4 -v

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


@timer("Finish {} in".format(__name__.rsplit(".", maxsplit=1)[0]))
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

    # Pretend to empty mailbox
    pretend_to_empty_mailbox()

    # Start an interactive session
    if util.check_options("-I", "--interactive"):
        from where.tools import interactive

        interactive.interactive(rundate, pipeline)
        return

    # Set up a new analysis or update an existing one
    setup.setup(rundate, pipeline)

    # Run the analysis
    with timer(f"Finish pipeline {pipeline.upper()} in"):
        pipelines.run(rundate, pipeline)


def pretend_to_empty_mailbox():
    """Pretend to try to empty a mailbox

    This function simply prints an error message to stderr. No actual attempt to empty a mailbox is done. This is
    included purely for nostalgic (and possibly some misunderstood backwards compatibility) reasons :P
    """
    if util.check_options("-M", "--mailbox"):
        print("/var/spool/mail/pha: Permission denied.", file=sys.stderr)


def profiler():
    """Run Where with profiling turned on

See below for regular options for running Where.

Added options when profiling:

--line_profile=f1,...   Do line profiling of functions <f1>, <...>.
--show_profile          Print profile information to the terminal after running.
--show_profile=col:num  Show <num> entries sorted on column <col> when printing.
--profile_output=file   Store profile output in <file>.

If no profiling options are given, the default behaviour is do profiling on the
function level and to not print anything to the terminal. The default filenames
are `where.lprof` for line profiling and `where.prof` for function profiling.
The profile output will be stored in the current directory.


line_profiler:
--------------

When doing line profiling, the functions to profile should be specfied with
their full name (possibly without the `where.`-prefix), for instance

    {exe:profiler} 2009 11 2 --line_profile=where.__main__.parse_techs

Several functions can be line profiled by separating them by comma,

    {exe:profiler} 2015 8 4 --line_profile=__main__.parse_techs,lib.log.init

It is also possible to do simple wildcard matching with a * (only at the end of
the function name). For instance will the following profile all functions in
`lib.log` starting with `l`:

    {exe:profiler} 2016 3 1 --line_profile=lib.log.l*

The results from doing a line profile can be printed using:

    python -m line_profiler where.lprof

Note that the data file does not store the actual source code, so if the source
file has been changed since the profiling was run, the print out will most
likely be wrong.


KCacheGrind:
------------

To analyse regular profiling data, one possibility is to use KCacheGrind
(kcachegrind.github.io). A simple way to start KCacheGrind is with the
pyprof2calltree-utility:

    pyprof2calltree -k -i where.prof

If these utilities are not installed use

    sudo apt-get install kcachegrind
    conda install pyprof2calltree

to install them.

See https://docs.kde.org/trunk5/en/kdesdk/kcachegrind/kcachegrind.pdf for more
information on KCacheGrind.


pstats:
-------

A simpler option is to use the pstats-module in the Python standard library. See
https://docs.python.org/3/library/profile.html for the documentation.

One interesting option is the interactive Profile Statistics Browser that can be
started by:

    python -m pstats where.prof

Type 'help' to see a list of commands, and 'help <command>' for information on a
particular command. Some simple examples are

    stats 10                        # Show statistics for 10 entries
    strip                           # Shorter filenames
    sort time                       # Sort on internal time
    callers unique                  # Show callers of functions matching unique
    callees azimuth                 # Show which functions azimuth calls


Where:
------
    """
    # Use profiler doc for help message when running where_profiler
    sys.modules[__name__].__doc__ = "{}\n\n{}".format(
        sys.modules[__name__].profiler.__doc__.strip(),
        sys.modules[__name__].__doc__.replace("{exe}", "{exe:profiler}"),
    )

    # Set up profiler and call main as usual
    from where.lib import profiler

    profiler.profiler(main)


# Run main function only when running as script
if __name__ == "__main__":
    sys.exit(main())
