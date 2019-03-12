"""Where library module for logging

Description:
------------

This module provides simple logging inside Where. To write a log message, simply call one of where.log-functions
corresponding to the log levels defined in where.lib.enums.


Example:
--------

    >>> from where.lib import log
    >>> log.init("info", prefix="My prefix")
    >>> n, m = 5, 3
    >>> log.info(f"Calculating the inverse of a {n:>2d}x{m:<2d} matrix")
    INFO  [My prefix] Calculating the inverse of a  5x3  matrix

"""

# Standard library imports
import functools

# Midgard imports
from midgard.dev import log as mg_log

# Where imports
from where.lib import enums  # Log levels and colors for Where
from where.lib import exceptions

# Make functions from Midgard available
from midgard.dev.log import log, blank, init, file_init, print_file  # noqa


# Make each log level available as a function, done here to include extra Where log levels
for level in enums.get_enum("log_level"):
    globals()[level.name] = functools.partial(mg_log.log, level=level.name)


# Overwrite log.fatal to raise an exception
def fatal(log_text):
    mg_log.log(log_text, "fatal")
    raise exceptions.WhereExit(f"Exiting Where due to {log_text!r}") from None
