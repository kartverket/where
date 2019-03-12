"""Run Where with profiling turned on

Usage:

    {exe:profiler} [arguments]

Typically, arguments and options are the same as when running Where
normally. See below for regular options for running Where.

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

# Standard library imports
import sys

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import __main__
from where.lib import profiler as _profiler


@plugins.register
def profiler():
    """Run Where with profiling turned on"""
    # Use profiler doc for help message when running where_profiler
    sys.modules[__main__.__name__].__doc__ = "{}\n\n{}".format(
        sys.modules[__name__].__doc__.strip(),
        sys.modules[__main__.__name__].__doc__.replace("{exe}", "{exe:profiler}"),
    )

    # Set up profiler and call main as usual
    _profiler.profiler(__main__.main)
