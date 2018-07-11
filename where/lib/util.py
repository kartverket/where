"""Where library module with utility functions for easier script development

Example:
from where.lib import util
directory, date = util.parse_args('string', 'date')

Description:

This module provides the boilerplate code necessary for starting a
script. In particular handling of command line arguments and default
options including --help are done.

"""

# Standard library imports
from datetime import date, datetime
import functools
import getpass
import os.path
import platform
import re
import sys

# Where imports
import where
from where.lib import config
from where.lib import files
from where.lib import log

# Source code cache, used by trace_source
_CACHE_SRC = dict()


def read_option_value(option, default=""):
    """Read the value of one command line option

    The option should be specified as a string with the necessary - or -- in front. If that option is not one of the
    command line arguments, default is returned. If there is a value following the option that value is returned as a
    string (separated by =). If there are several occurences of the option, the first one is returned.

    Args:
        option (String):    Option specified with the leading - or --.
        default (String):  Optional default value that is returned if the option is not specified.

    Returns:
        String: The option or the value of the option. The default value if the option is not specified.
    """
    for arg in sys.argv[1:]:
        if arg.startswith(f"{option}="):
            # Remove the part up to and including the first =-sign
            return arg.split("=", maxsplit=1)[-1]

    return default


def check_options(*options):
    """Check if any of a list of options is specified on the command line

    Returns the actual option that is specified. The first option specified on the command line is returned if there
    are several matches. Returns the empty string if no option is specified. This means that this method works fine
    also in a boolean context, for example

        if check_options('-F', '--force'):
            do_something()

    Args:
        options:   Strings specifying which options to check for, including '-'-prefix.

    Returns:
        String: Option that is specified, blank string if no option is specified

    """
    cmd_argv = [a.split("=")[0] for a in sys.argv[1:]]

    for option in cmd_argv:
        if option in options:
            return option

    return ""


def no_traceback(func):
    """Decorator for turning off traceback, instead printing a simple error message

    Use the option --show_tb to show the traceback anyway.
    """
    if check_options("-T", "--showtb"):
        return func

    def no_traceback_hook(_not_used_1, value, _not_used_2):
        """Only prints the error message, no traceback."""
        log.error(str(value))

    @functools.wraps(func)
    def _no_traceback(*args, **kwargs):
        remember_excepthook = sys.excepthook
        sys.excepthook = no_traceback_hook
        values = func(*args, **kwargs)
        sys.excepthook = remember_excepthook
        return values

    return _no_traceback


@no_traceback
def parse_args(*param_types, doc_module=None):
    """Parse command line arguments and general options

    Show the help message parsed from the script's docstring if -h or
    --help option is given. Show the script's version if --version is
    given. Log versions of python, the script and the configuration.
    Finally parse arguments from the given parameter types.

    Args:
        param_types: Strings describing the expected parameter types.
                     Each string must be one of the keys in #_PARSERS.

    Returns:
        List of command line arguments parsed according to param_types.
    """
    # Help message
    if check_options("-h", "--help"):
        _print_help_from_doc(doc_module)
        sys.exit(0)

    # Version information
    if check_options("--version"):
        _print_version_from_doc(doc_module)
        sys.exit(0)

    # Log version of python and the program, and the configuration file used
    if doc_module:
        log.info("Start {} at {}", _get_program_version(doc_module), datetime.now().strftime(config.FMT_datetime))
        log.python_version()
        log.debug("Receive command line arguments [{}]", ", ".join(sys.argv[1:]))
        log.configuration(cfg=get_program_name())

    # Parse arguments
    try:
        arguments = [_PARSERS[type]() for type in param_types]
    except Exception:
        _print_help_from_doc(doc_module)
        raise

    # Return arguments (scalar if only one element, None if list is empty)
    if len(arguments) > 1:
        return arguments
    elif arguments:
        return arguments[0]


def not_implemented():
    """A placeholder for functions that are not implemented yet

    A note about the missing implementation is written to the log.
    """
    caller = sys._getframe(1)
    funcname = caller.f_code.co_name
    args = ", ".join([k + "=" + str(v) for k, v in caller.f_locals.items()])
    filename = caller.f_code.co_filename
    lineno = caller.f_lineno
    log.fatal("Function {}({}) is not implemented in {}, line {}", funcname, args, filename, lineno)


def get_callers():
    """Get a list of methods calling this function

    """
    callers = list()
    caller = sys._getframe(1)
    while caller:
        func_name = caller.f_code.co_name
        line_no = caller.f_lineno
        file_name = caller.f_code.co_filename
        if file_name.endswith(".py"):
            module_name = file_name[-file_name[::-1].find("/erehw/"):-3].replace("/", ".")
            callers.insert(0, "{}.{} ({})".format(module_name, func_name, line_no))

        caller = caller.f_back

    return "\n    -> ".join(callers)


def get_program_name():
    """Get the name of the running program

    Returns:
        String trying to be similar to how the user called the program.
    """
    program_name = sys.argv[0]
    if not program_name.startswith("./"):
        program_name = os.path.basename(program_name)
    return program_name


def get_program_info():
    """Get the name and the version of the running program
    """
    return f"{get_program_name()} v{where.__version__}"


def get_authors_doc(*authors):
    """Get list of authors of software

    Returns:
        String:  Name and contact for authors.
    """
    return "\n".join(
        "* {author.name}, <{author.email}>".format(author=a) for a in where._AUTHORS if a.start < date.today() < a.end
    )


def get_user_info(user=None):
    """Get info about the user running the software

    The info must be registered in the configuration file in the `user_info`-section.

    Args:
        user (String):  Optional user name, if not given the current user is used

    Returns:
        Dict:  Information about user
    """
    user = config.program.get("user", value=user, default=getpass.getuser()).str
    user_info = config.where.get(section="user_info", key=user, default="").as_tuple(", *")

    info_dict = dict(user=user, **dict(zip(("name", "email", "inst_abbreviation"), user_info)))
    if "inst_abbreviation" in info_dict:
        inst_info = config.where.get(
            section="institution_info", key=info_dict["inst_abbreviation"].lower(), default=""
        ).as_tuple(
            ", *"
        )
        info_dict.update(dict(zip(("inst_name", "inst_address"), inst_info)))

    return info_dict


def write_requirements():
    """Write requirements (python modules) to file for reproducibility.

    Note that this only stores the modules that have been imported, and that have a `__version__`-attribute (see PEP
    396 - https://www.python.org/dev/peps/pep-0396/)
    """
    # Find versions of imported modules (use list() for copy in case modules are imported when reading __version__)
    reqs = {n: getattr(m, "__version__", None) for n, m in list(sys.modules.items())}
    reqs["python"] = platform.python_version()
    reqs_str = "\n".join(sorted("{}=={}".format(m, v.strip()) for m, v in reqs.items() if isinstance(v, str)))

    # Write to requirements-file
    with files.open("requirements", mode="w") as fid:
        fid.write(reqs_str + "\n")


def trace_full(frame, event, arg):
    """Show full trace for a given stack frame

    This function is mainly meant for debug and can for instance be activated by inserting either
    `sys.settrace(util.trace_full)` or `sys.setprofile(util.tracefull)` in the source code where tracing should start.

    Args:
        frame:    Current frame object.
        event:    Event, not used but needed to be compatible with sys.setprofile.
        arg:      Arg, not used but needed to be compatible with sys-.setprofile.
    """
    callers = list()
    while frame:
        func_name = frame.f_code.co_name
        line_no = frame.f_lineno
        file_name = frame.f_code.co_filename
        if file_name.endswith(".py"):
            module_name = file_name[-file_name[::-1].find("/erehw/"):-3].replace("/", ".")
            callers.insert(0, "{}.{} ({})".format(module_name, func_name, line_no))

        frame = frame.f_back

    print("\n    -> ".join(callers))


def trace_source(frame, event, arg):
    """Show full trace for a given stack frame

    This function is mainly meant for debug and can for instance be activated by inserting
    `sys.setprofile(util.tracefull)` in the source code where tracing should start. This function will then be called
    whenever a function is called or returned from. It is also possible to use `sys.settrace` in a similar way.

    Args:
        frame:             Current frame object.
        event (String):    Event, not used but needed to be compatible with sys.setprofile.
        arg:               Arg, not used but needed to be compatible with sys-.setprofile.
    """
    caller = frame.f_back
    if caller is None:
        return

    file_name = caller.f_code.co_filename
    if "/where/" not in file_name:
        return

    if file_name not in _CACHE_SRC:
        with open(file_name, mode="r") as fid:
            _CACHE_SRC[file_name] = {no: ln.strip() for no, ln in enumerate(fid.readlines(), start=1)}

    line_no = caller.f_lineno
    module_name = file_name[-file_name[::-1].find("/erehw/"):-3].replace("/", ".")
    func_name = "{}.{} ({}):".format(module_name, caller.f_code.co_name, line_no)

    print("-> {:<40s} {}".format(func_name, _CACHE_SRC[file_name][line_no]))


def _get_doc(doc_module=None):
    """Get the docstring of the running program

    Args:
        doc_module:  String, name of the module with docstring. Default is __main__.

    Returns:
        String, the docstring of the given module.
    """
    if doc_module is None:
        doc_module = "__main__"

    doc = sys.modules[doc_module].__doc__

    return "" if doc is None else doc


def _get_program_version(doc_module=None):
    """Get version info for the running program from SVN keywords

    The version of the program is read from the main __init__.py-file of Where. In addition, the docstring of the
    running script is parsed for revision information. In particular we look for the SVN keyword substitutions Revision
    and Date. If no revision information is found, only the program name and version is returned.

    Args:
        doc_module:  String, name of the module to be searched for revision information. Default is __main__.

    Returns:
        String, information about the version of the given module (running script).
    """
    doc = _get_doc(doc_module)
    revision = re.search(r"\$Revision: (\d+) \$", doc)
    last_changed = re.search(r"\$Date: (\d{4}-\d{2}-\d{2}) ", doc)
    if revision and last_changed:
        revision_text = ", revision {} ({})".format(revision.groups()[0], last_changed.groups()[0])
    else:
        revision_text = ""
    return "{} v{}{}".format(get_program_name(), where.__version__, revision_text)


def _print_help_from_doc(doc_module=None):
    """Filter the docstring to make it a better help text and print it

    Removes @-directives in the docstring. Furthermore, we replace the
    Example-heading used by Doxygen with Usage as that makes more
    sense in a help text. Finally $-symbols at the beginning and end
    of lines, as introduced by the SVN keyword substitution, are
    removed.
    """
    from where import pipelines  # Local import to avoid circular import

    replace_vars = dict(
        pipelines_doc=pipelines.doc(),
        maintainers="{maintainers}",  # Handled by where._update_doc
        version=where.__version__,
        exe=where.__executable__,
    )
    doc = where._update_doc(_get_doc(doc_module).format(**replace_vars))

    for line in doc.splitlines():
        if line.startswith("@"):
            continue
        line = re.sub(r"::", ":", line)
        line = re.sub(r"^\$| ?\$$", "", line)
        print(line)


def _print_version_from_doc(doc_module=None):
    """Use docstring information to print a version string

    Uses the $-keywords from SVN to create a string containing version
    information.
    """
    print(_get_program_version(doc_module))


def _next_argument():
    """Return the next argument as a string

    The next argument is returned. Options (strings starting with -)
    are ignored. The argument is removed from the sys.argv-list.

    Returns:
        String with the next argument.
    """
    for idx in range(1, len(sys.argv)):
        if not sys.argv[idx].startswith("-"):
            return sys.argv.pop(idx)

    raise TypeError("Missing argument")


def _parse_date():
    """Return the three next arguments as a date

    The three next arguments are parsed as a datetime-object. An
    error-message is raised if the arguments are not of the required
    form.

    Returns:
        Datetime-object with the date specified in the three next arguments.
    """
    try:
        return date(_parse_int(), _parse_int(), _parse_int())
    except (ValueError, TypeError) as err:
        err.args = (f"{err.args[0]}\n  Date should be written year month day, for instance {date.today():%Y %-m %-d}",)
        raise


def _parse_doy():
    """Return the two next arguments as a date

    The two next arguments are parsed as a year and day-of-year. An error-message is raised if the arguments are not of
    the required form.

    Returns:
        Datetime: The date specified as year and day-of-year in the two next arguments.
    """
    try:
        return datetime.strptime("{} {}".format(_parse_int(), _parse_int()), "%Y %j").date()
    except ValueError as err:
        err.args = (err.args[0] + "\n    Date should be written year doy, for instance 2009 306",)
        raise


def _parse_float():
    """Return the next argument as a float

    The next argument is returned. Options (strings starting with -)
    are ignored. The argument is removed from the sys.argv-list. A
    ValueError is raised if the argument can not be parsed as a float.

    Returns:
        Float with the next argument.
    """
    return float(_next_argument())


def _parse_int():
    """Return the next argument as an int

    The next argument is returned. Options (strings starting with -)
    are ignored. The argument is removed from the sys.argv-list. A
    ValueError is raised if the argument can not be parsed as an int.

    Returns:
        Int with the next argument.
    """
    return int(_next_argument())


# Parsers for each parameter type
_PARSERS = {"date": _parse_date, "doy": _parse_doy, "float": _parse_float, "int": _parse_int, "string": _next_argument}

# Store string about how program is started
COMMAND = "{} {}".format(get_program_name(), " ".join(sys.argv[1:]))
