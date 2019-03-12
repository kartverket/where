"""Run extra tools related to Where

Usage:

    {exe:tools} <tool> [arguments]

The program requires the name of a tool. Each tool has its own set of
arguments. Run "{exe:tools} <tool> -h" to get more information and help about a
particular tool.

The following tools are available:

===================  ===========================================================
Tool                 Description
===================  ===========================================================
{tools_doc}
===================  ===========================================================

The following general options are recognized for all tools:

===================  ===========================================================
Option               Description
===================  ===========================================================
--debug, ...         Show additional debug information. Other flags such as
                     --all, --debug, --time, --dev, --info, --out, --warn,
                     --check, --error, --fatal, --none are also allowed, and
                     shows differing amounts of information as the program runs.
--version            Show version information and exit.
-h, --help           Show this help message and exit.
-T, --showtb         Show traceback if the tool crashes.
===================  ===========================================================


Current Maintainers:
--------------------

{maintainers}

Version: {version}

"""

# Standard library imports
import sys

# Midgard imports
from midgard.dev import exceptions
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log
from where.lib import util
from where import pipelines


@util.no_traceback
def main():
    """Invoke where_tools

    To add a new tool, simply add a new .py-file with a registered plugin that
    will be called when the tool is called.

    """
    # Start logging
    log.init(log_level="info")

    # First read tool from command line, don't use util.parse_args() as that overrides -h for each tool
    try:
        tool = [a for a in sys.argv[1:] if not a.startswith("-")][0]
        sys.argv.remove(tool)
    except IndexError:
        util._print_help_from_doc(__name__)
        raise SystemExit

    # Check that tool is available and figure out signature
    try:
        sig = plugins.signature(__name__, tool)
        tool_module = plugins.get(__name__, tool).function.__module__
    except exceptions.UnknownPluginError as err:
        util._print_help_from_doc(__name__)
        err.args = (f"{err.args[0]}\n    Available tools are {', '.join(plugins.names(__name__))}",)
        raise

    # Parse parameters
    util.check_help_and_version(doc_module=tool_module)

    tool_args = dict()
    for key, param in sig.parameters.items():
        if param.annotation is None:
            raise SystemExit(f"{param} in {tool} tool is not annotated")

        if param.annotation == "datedoy":
            if util.check_options("--doy"):
                date = util.parse_args("doy")  # TODO: Does not work ("doy", doc_module=__name__)
            else:
                date = util.parse_args("date")  # TODO: Does not work ("date", doc_module=__name__)
            tool_args[key] = date

        elif param.annotation == "pipeline":
            tool_args[key] = pipelines.get_from_options()
            config.read_pipeline(tool_args[key])

        elif param.annotation == "option":
            tool_args[key] = util.read_option_value(f"--{key}")

        else:
            tool_args[key] = util.parse_args(param.annotation, doc_module=tool_module)

    # Call tool
    plugins.call(__name__, tool, **tool_args)


def doc():
    return "\n".join(f"{t:<20} {d}" for t, d in plugins.doc_all(__name__, long_doc=False, use_module=True).items())
