"""Where library module for logging

Example:
--------

from where.lib import log
log.init()
log.info('Calculating the inverse of a {:>2d}x{:<2d} matrix', n, m)

Description:
------------

This module provides simple logging inside Where. To write a log message, simply call one of where.log.debug,
where.log.info, where.log.warn, where.log.error or where.log.fatal with a log message that will be parsed through
str.format.

In the future we may look into wrapping the standard library logging module instead of writing our own.
"""

# Standard library imports
import atexit
from datetime import datetime
import os.path
import platform
import shutil
import sys

# Where imports
from midgard.dev import console

from where.lib import config
from where.lib import exceptions

# Log-levels, a log message will be shown if the message log level is at or above the current log level.
_LEVELS = ("ALL", "DEBUG", "TIME", "DEV", "INFO", "OUT", "WARN", "CHECK", "ERROR", "FATAL", "NONE")
LOGLEVELS = {level: num for num, level in enumerate(_LEVELS)}

# Colors used when printing log, log-levels without an entry will not be colored.
LOGCOLORS = {
    "DEV": console.color.Fore.BLUE,
    "TIME": console.color.Fore.WHITE,
    "OUT": console.color.Style.BRIGHT,
    "CHECK": console.color.Style.BRIGHT + console.color.Fore.YELLOW,
    "WARN": console.color.Fore.YELLOW,
    "ERROR": console.color.Fore.RED,
    "FATAL": console.color.Style.BRIGHT + console.color.Fore.RED,
}

# Information about the state of logging. Updated in the init- and file_init-functions.
LOGINFO = {"level": None, "do_console": False, "do_file": False, "keep_cache": True}

# Before logging is started explicitly, keep a cache of logs in memory
LOGCACHE = list()


def log(
    level,
    log_text,
    *fmtargs,
    log_to_console=True,
    log_to_file=True,
    cache_data=None,
    keep_cache=True,
    print_level=True,
    print_caller=False,
    **fmtkws,
):
    """Write message to log on both console and file

    The message is parsed with str.format and written to the log if the message log level id at or above the current
    log level.  Typically you should call one of log.debug, log.info, log.warn, log.error or log.fatal instead of
    calling this function directly.

    Args:
        level (String):           Message log level, one of the keys of LOGLEVELS.
        log_text (String):        String, will be parsed with str.format.
        fmtargs (List):           Positional arguments passed on to str.format.
        log_to_console (Boolean): Write message to console.
        log_to_file (Boolean):    Write message to file.
        cache_data (Dict):        Data from cache. Only used by _dump_cache().
        keep_cache (Boolean):     Store log message in cache if it is not logged to console and file.
        print_level (Boolean):    Prefix log message with log level.
        print_caller (Boolean):   Print information about calling function (Always done if log-level is DEBUG).
        fmtkws (Dict):            Keyword arguments passed on to str.format.
    """
    timestamp = datetime.now().strftime(config.FMT_datetime) if cache_data is None else cache_data["timestamp"]
    text = log_text.format(*fmtargs, **fmtkws)
    caller = sys._getframe(2) if cache_data is None else cache_data["caller"]
    try:
        runinfo = "[{} {} {}] ".format(
            config.analysis.analysis.str.upper(), config.analysis.session.str, config.analysis.rundate.str
        )
    except exceptions.MidgardException:
        runinfo = ""

    # Add more information if we are logging at the DEBUG level
    if LOGINFO["level"] == LOGLEVELS["DEBUG"] or print_caller:
        func_name = caller.f_code.co_name
        file_name = caller.f_code.co_filename
        line_num = caller.f_lineno
        text = "{}\n         [{} ({}:{})]".format(text, func_name, file_name, line_num, text)

    # Store log in cache if we are not already logging to both console and file
    if keep_cache and LOGINFO["keep_cache"] and not (LOGINFO["do_console"] and LOGINFO["do_file"]):
        LOGCACHE.append(
            dict(
                timestamp=timestamp,
                level=level,
                text=text,
                caller=caller,
                to_console=log_to_console & (not LOGINFO["do_console"]),
                to_file=log_to_file & (not LOGINFO["do_file"]),
            )
        )

    # First do the console-log
    if log_to_console and LOGINFO["do_console"] and LOGLEVELS.get(level, max(LOGLEVELS.values())) >= LOGINFO["level"]:
        is_warning = LOGLEVELS.get(level, max(LOGLEVELS.values())) >= LOGLEVELS["WARN"]
        stream = sys.stderr if is_warning else sys.stdout
        color = LOGCOLORS.get(level, "")
        level_str = f"{color}{level + ':':<6}" if print_level else ""
        print(f"{level_str} {runinfo}{text}", file=stream)

    # Next do the file log
    if log_to_file and LOGINFO["do_file"] and LOGLEVELS.get(level, max(LOGLEVELS.values())) >= LOGINFO["file_level"]:
        level_str = f"{level + ':':<6}" if print_level else ""
        LOGINFO["file_id"].write(f"{level_str} {runinfo}{text}\n")


def blank(log_to_console=True, log_to_file=False, keep_cache=True):
    """Write a blank line to the console log to make it easier to read

    The file log is primarily meant to be machine-readable, so by default we do not log the blank line to file. The
    blank line is written as if at the INFO level.

    Args:
        log_to_console: Write message to console.
        log_to_file:    Write message to file.
    """
    if keep_cache and not (LOGINFO["do_console"] and LOGINFO["do_file"]):
        LOGCACHE.append(
            dict(
                timestamp=None,
                level=None,
                text=None,
                to_console=log_to_console & (not LOGINFO["do_console"]),
                to_file=log_to_file & (not LOGINFO["do_file"]),
            )
        )

    if log_to_console and LOGINFO["do_console"] and LOGLEVELS["INFO"] >= LOGINFO["level"]:
        print("", file=sys.stdout)

    if log_to_file and LOGINFO["do_file"] and LOGLEVELS["INFO"] >= LOGINFO["file_level"]:
        LOGINFO["file_id"].write("\n")


def debug(log_text, *fmtargs, **fmtkws):
    """Write debug message to log

    The message is parsed with str.format and written to the log if the current log level is DEBUG or lower.

    Args:
        log_text:   String, will be parsed with str.format.
        fmtargs:    Arguments passed on to str.format.
        fmtkws:     Keyword arguments passed on to str.format.
    """
    log("DEBUG", log_text, *fmtargs, **fmtkws)


def murks(log_text, *fmtargs, **fmtkws):
    """Write temporary debug message to log

    The message is parsed with str.format and only written to the log if `--murks` is specified on the command
    line. The message is not written to the file log.

    It is also possible to specify `--murks=prefix` on the command line. In this case, if the message starts with
    `prefix` an interactive session will be started after the murks message is logged.

    Args:
        log_text:   String, will be parsed with str.format.
        fmtargs:    Arguments passed on to str.format.
        fmtkws:     Keyword arguments passed on to str.format.
    """
    from where.lib import util

    if not util.check_options("--murks"):
        return

    # Log murks messages
    log("MURKS", log_text, *fmtargs, log_to_file=False, **fmtkws)

    # Start an interactive session (with variables from correct namespace) if message matches command line prefix
    interactive_prefix = util.read_option_value("--murks")
    if interactive_prefix is not None and log_text.startswith(interactive_prefix):
        import IPython

        frame = sys._getframe(1)
        namespace = dict(**frame.f_globals, **frame.f_locals)
        IPython.embed(user_ns=namespace)


def dev(log_text, *fmtargs, **fmtkws):
    """Write dev message to log

    The message is parsed with str.format and written to the log if the current log level is DEV or lower.

    Args:
        log_text:   String, will be parsed with str.format.
        fmtargs:    Arguments passed on to str.format.
        fmtkws:     Keyword arguments passed on to str.format.
    """
    log("DEV", log_text, *fmtargs, **fmtkws)


def time(log_text, *fmtargs, **fmtkws):
    """Write time message to log

    The message is parsed with str.format and written to the log if the current log level is TIME or lower.

    Args:
        log_text:   String, will be parsed with str.format.
        fmtargs:    Arguments passed on to str.format.
        fmtkws:     Keyword arguments passed on to str.format.
    """
    log("TIME", log_text, *fmtargs, **fmtkws)


def info(log_text, *fmtargs, **fmtkws):
    """Write info message to log

    The message is parsed with str.format and written to the log if the current log level is INFO or lower.

    Args:
        log_text:   String, will be parsed with str.format.
        fmtargs:    Arguments passed on to str.format.
        fmtkws:     Keyword arguments passed on to str.format.
    """
    log("INFO", log_text, *fmtargs, **fmtkws)


def out(log_text, *fmtargs, **fmtkws):
    """Write output message to log

    The message is parsed with str.format and written to the log if the current log level is INFO or lower.

    Args:
        log_text:   String, will be parsed with str.format.
        fmtargs:    Arguments passed on to str.format.
        fmtkws:     Keyword arguments passed on to str.format.
    """
    log("OUT", log_text, *fmtargs, **fmtkws)


def check(log_text, *fmtargs, **fmtkws):
    """Write message to log about something that should be checked

    The message is parsed with str.format and written to the log if the current log level is CHECK or lower.

    Args:
        log_text:   String, will be parsed with str.format.
        fmtargs:    Arguments passed on to str.format.
        fmtkws:     Keyword arguments passed on to str.format.
    """
    log("CHECK", log_text, *fmtargs, **fmtkws)


def warn(log_text, *fmtargs, **fmtkws):
    """Write warning message to log on stderr

    The message is parsed with str.format and written to the log if the current log level is WARN or lower.

    Args:
        log_text:   String, will be parsed with str.format.
        fmtargs:    Arguments passed on to str.format.
        fmtkws:     Keyword arguments passed on to str.format.
    """
    log("WARN", log_text, *fmtargs, **fmtkws)


def error(log_text, *fmtargs, **fmtkws):
    """Write error message to log on stderr

    The message is parsed with str.format and written to the log if the current log level is ERROR or lower.

    Args:
        log_text:   String, will be parsed with str.format.
        fmtargs:    Arguments passed on to str.format.
        fmtkws:     Keyword arguments passed on to str.format.
    """
    log("ERROR", log_text, *fmtargs, **fmtkws)


def fatal(log_text, *fmtargs, print_caller=True, **fmtkws):
    """Write fatal message to log on stderr

    The message is parsed with str.format and written to the log if the current log level is FATAL or lower. The
    message is not written to console using normal log functionality. Instead, a SystemExit exception is raised with
    the given message. This also ends the running script.

    Args:
        log_text (String):        Log-text, will be parsed with str.format.
        fmtargs:                  Arguments passed on to str.format.
        print_caller (Boolean):   Print information about calling function.
        fmtkws:                   Keyword arguments passed on to str.format.
    """
    log("FATAL", log_text, print_caller=print_caller, *fmtargs, **fmtkws)
    raise exceptions.WhereExit("Exiting Where due to '{}'.".format(log_text.format(*fmtargs, **fmtkws))) from None


def assert_true(assertion, log_text, *fmtargs, **fmtkws):
    """Test that an assertion is True, exit the program if it is not.

    Args:
        assertion (Bool):    Assertion, the program will exit if the assertion is False.
        log_text (String):   Log-text, will be parsed with str.format.
        fmtargs:             Arguments passed on to str.format.
        fmtkws:              Keyword arguments passed on to str.format.
    """
    if not assertion:
        fatal(log_text, *fmtargs, print_caller=False, **fmtkws)


def assert_not_none(obj, log_text, *fmtargs, **fmtkws):
    """Test that an object is not None, exit the program if it is.

    Args:
        obj (Object):        Object to test, the program will exit if the object is None.
        log_text (String):   Log-text, will be parsed with str.format.
        fmtargs:             Arguments passed on to str.format.
        fmtkws:              Keyword arguments passed on to str.format.
    """
    if obj is None:
        fatal(log_text, *fmtargs, print_caller=False, **fmtkws)


def lowest(*loggers):
    """Choose the logger at the lowest level

    Args:
        loggers (List):  Log functions like debug(), info(), etc.

    Returns:
        Function:   The log function with lowest log level.
    """
    log_levels = {LOGLEVELS[f.__name__.upper()]: f for f in loggers}
    return min(log_levels.items())[1]


def highest(*loggers):
    """Choose the logger at the highest level

    Args:
        loggers (List):  Log functions like debug(), info(), etc.

    Returns:
        Function:   The log function with highest log level.
    """
    log_levels = {LOGLEVELS[f.__name__.upper()]: f for f in loggers}
    return max(log_levels.items())[1]


def python_version():
    """Log info message with information about the python interpreter

    Use sys.executable and sys.version_info to find name and version of the python interpreter. Use the platform.uname
    to find servername instead of os.uname because the latter is not supported on Windows.
    """
    pybin = os.path.basename(sys.executable)
    debug(
        "Use {pybin}, version {v[0]}.{v[1]}.{v[2]} on process {pid}@{server}",
        pybin=pybin,
        v=sys.version_info[0:3],
        pid=os.getpid(),
        server=platform.uname().node,
    )


def configuration(cfg="where"):
    """Log info message with information about the configuration file

    Finds the name by using the where.config module.

    Args:
        cfg:  Which configuration to report files from.
    """
    try:
        config_files = getattr(config, cfg).sources
    except AttributeError:
        pass
    else:
        info(f"Use {cfg.title()} configuration from {', '.join(config_files)}")


def init(log_level=None):
    """Parse command line options and set the current log level accordingly

    Look through the command line arguments and set the current log level to the first (if any) option of the form
    --level where level is one of the keys in LOGLEVELS (case-insensitive). If no command line arguments are found,
    the log:default_level value from the config value is used. If this does not exist, the log level is set to INFO.

    Args:
        log_level (String):  Optional log level, default is to find the log level in config or command line arguments.
    """
    LOGINFO["do_console"] = True

    # Default choice: set log level to INFO
    log_levels = ["INFO"]

    # 3rd priority: log level from config file
    log_levels.append(config.where.log.default_level.str.upper())

    # 2nd priority: log level from command line argument
    for arg in sys.argv:
        if arg.startswith("--"):
            log_levels.append(arg[2:].upper())

    # 1st priority: log level from function argument
    if log_level is not None:
        log_levels.append(log_level.upper())

    # Set log level to last valid choice
    LOGINFO["level_name"] = [level for level in log_levels if level in LOGLEVELS].pop()
    if LOGLEVELS[LOGINFO["level_name"]] != LOGINFO["level"]:
        LOGINFO["level"] = LOGLEVELS[LOGINFO["level_name"]]
        info("Set log-level to {}", LOGINFO["level_name"])

    _dump_cache()


def file_init(log_path=None, log_level=None, append=False):
    """Set up and start the logging to file

    Reads the Where configuration and starts the file logging if log:log_to_file is True. The proper log level is also
    set from the configuration. The log file is set up so that it will close automatically when the script ends.

    Args:
        log_path (String):  Optional path to log file. Default is to find log path in config.
        log_level (String): Optional log level, default is to find the log level in config or command line arguments.
        append (Boolean):   Append to log file or start new (with log rolling).
    """
    if log_path is None:
        if not config.where.log.log_to_file.bool:
            return

    # Close any file currently being logged
    if LOGINFO["do_file"]:
        atexit.unregister(file_end)
        file_end()

    # Set the logging level
    log_levels = [
        "TIME",
        LOGINFO["level_name"],
        config.where.log.filelog_level.str.upper(),
        None if log_level is None else log_level.upper(),
    ]
    LOGINFO["file_level"] = LOGLEVELS[[level for level in log_levels if level in LOGLEVELS].pop()]

    # Open the log file
    LOGINFO["do_file"] = True
    if log_path is None:
        program = os.path.splitext(os.path.basename(sys.argv[0]))[0]
        log_directory = config.where.log.path.path
        log_path = log_directory / "{}.log".format(program)
    _setup_logfiles(log_path, do_log_rolling=not append)
    mode = "a" if append else "w"
    LOGINFO["file_id"] = open(log_path, mode=mode)
    atexit.register(file_end)  # Close log file when script exits

    # Send cache to file
    _dump_cache()


def file_end():
    """Close the logging to file
    """
    if not LOGINFO["do_file"]:
        return

    LOGINFO["do_file"] = False
    LOGINFO["file_id"].close()


def print_file(file_path, log_level="INFO"):
    if log_level not in LOGLEVELS:
        error("Log level '{}' is unknown. Use one of {}", log_level, ", ".join(LOGLEVELS))
        return

    info("Showing contents of log file {} at level {} and above", file_path, log_level)
    level_limit = LOGLEVELS[log_level]
    current_level = log_level
    with open(file_path, mode="r") as fid:  # lib.files is not available in lib.log
        for line in fid:
            line_level = line.split()[0].rstrip(":")
            current_level = line_level if line_level in LOGLEVELS else current_level
            if LOGLEVELS[current_level] >= level_limit:
                color = LOGCOLORS.get(current_level, "")
                print("{}| {}".format(color, line.rstrip()))


def turn_on_cache():
    """Turn on the caching of logs
    """
    LOGINFO["keep_cache"] = True


def turn_off_cache():
    """Turn off the caching of logs
    """
    LOGINFO["keep_cache"] = False
    LOGCACHE.clear()


def _dump_cache():
    """Dump the cache to console and/or file

    The cache is dumped to console and/or file. The cache is updated so that it will not be dumped to the same device
    twice.
    """
    for entry in LOGCACHE:
        if entry["level"] is None:
            blank(log_to_console=entry["to_console"], log_to_file=entry["to_file"], keep_cache=False)
        else:
            log(
                entry["level"],
                entry["text"],
                log_to_console=entry["to_console"],
                log_to_file=entry["to_file"],
                cache_data=dict(timestamp=entry["timestamp"], caller=entry["caller"]),
                keep_cache=False,
            )
        entry["to_console"] = False if LOGINFO["do_console"] else entry["to_console"]
        entry["to_file"] = False if LOGINFO["do_file"] else entry["to_file"]


def _setup_logfiles(log_path, do_log_rolling=True):
    """Perform necessary rolling of log files

    Creates the log directory if it does not already exist. Also rolls the log files. That is, if there are old log
    files, they will be moved to files with extension .0, .1 and so on. The number of rolled logs to keep is specified
    in the log:number_of_logs_backup field in the Where configuration.
    """
    log_directory = os.path.dirname(log_path)
    if not os.path.exists(log_directory):
        os.makedirs(log_directory)

    if not os.path.isfile(log_path):
        return

    num_rolls = config.where.log.number_of_log_backups.int if do_log_rolling else 0
    for roll in range(num_rolls - 1)[::-1]:
        if os.path.isfile("{}.{}".format(log_path, roll)):
            shutil.move("{}.{}".format(log_path, roll), "{}.{}".format(log_path, roll + 1))
    if num_rolls > 0:
        shutil.move(log_path, "{}.0".format(log_path))
