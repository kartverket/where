#!/usr/bin/env python3
"""Add a profiler when running Where





"""
# Standard library imports
import os

# Where imports
from where.lib import log
from where.lib import util


def profiler(func_to_profile):
    """Run a function with profiling turned on

    Args:
        func_to_profile (Function):   Function that will be called after profiling is turned on.
    """
    # Should we do line or function profiling?
    prof, info = _setup_line_profiler() if util.check_options("--line_profile") else _setup_cprofile()

    # Read command line options
    filename = util.read_option_value("--profile_output", default="where")
    if not os.path.splitext(filename)[1]:
        filename = f"{filename}.{info['extension']}"

    # Start profiling and run Where as usual. Store and print profile even if Where crashes
    log.info(f"Enable {info['doc']}, output stored in {filename}")
    prof.enable()
    try:
        func_to_profile()
    finally:
        prof.disable()

        # Store profile to file
        inspect_str = f"Inspect it using e.g. {' or '.join(info['inspect'])}"
        log.info(("Profile information stored to {f}\n       " + inspect_str).format(f=filename))
        prof.dump_stats(filename)

        # Print profile to terminal
        if util.check_options("--show_profile"):
            prof.print_stats()


def _setup_line_profiler():
    """Set up a line profiler"""
    import line_profiler

    prof = line_profiler.LineProfiler()
    funcs = util.read_option_value("--line_profile", default="").split(",")
    for fname in funcs:
        for func in _get_func_from_name(fname):
            prof.add_function(func)

    info = dict(
        extension="lprof",
        inspect=['"python -m line_profiler {f}"'],
        doc=f"line profiling functions {', '.join(funcs)}",
    )

    return prof, info


def _get_func_from_name(name):
    """Get a function object from a name

    Assume the name is given as module.function, or module.class.method.

    Args:
        name (String):   The full name of a function, for instance 'where.__main__._get_func_from_name'.

    Returns:
        Function:  The function object represented by name
    """
    import importlib
    import types

    # Import module and find function (possibly trying to find function inside class)
    name = name if name.startswith("where.") else "where." + name
    modname, _, fname = name.rpartition(".")
    try:
        obj = importlib.import_module(modname)
    except ImportError:
        modname, _, clsname = modname.rpartition(".")
        mod = importlib.import_module(modname)
        obj = getattr(mod, clsname)

    # Support simple *-wildcard matching
    if name.endswith("*"):
        fnames = [f for f in dir(obj) if f.startswith(fname[:-1]) and not f.startswith("__")]
        if not fnames:
            raise TypeError(f"Found no functions named {name!r}")
    else:
        fnames = [fname]

    # Generate each function in fnames
    for fname in fnames:
        func = getattr(obj, fname)

        # Unwrap decorated functions
        while hasattr(func, "__wrapped__"):
            func = func.__wrapped__

        # Return getter, setter and deleter functions of properties
        if isinstance(func, property):
            for accessor in ("fget", "fset", "fdel"):
                if getattr(func, accessor) is not None:
                    yield getattr(func, accessor)
        else:
            if isinstance(func, types.FunctionType):
                yield func
            elif not name.endswith("*"):
                raise TypeError(f"'{modname}.{fname}' is not a function")


def _setup_cprofile():
    """Set up a profiler"""
    import cProfile

    prof = cProfile.Profile()
    prof.print_stats = _print_cprofile_stats(prof)

    info = dict(extension="prof", inspect=['"pyprof2calltree -k -i {f}"', '"python -m pstats {f}"'], doc="profiling")

    return prof, info


def _print_cprofile_stats(prof):
    """Print statistics about data from cProfile

    Print statistics as defined by the `--show_profile` command line option. See the `profile`-function for more
    information.
    """
    import pstats

    profile_info = util.read_option_value("--show_profile", default="time:50")

    def print_stats():
        stat = pstats.Stats(prof).strip_dirs()
        if ":" in profile_info:
            sort_columns, _, amount = profile_info.partition(":")
        else:
            sort_keys = stat.get_sort_arg_defs().keys()
            sort_columns, amount = (profile_info, "") if profile_info in sort_keys else ("time", profile_info)

        amount = int(amount) if amount.isnumeric() else amount
        pstats.Stats(prof).strip_dirs().sort_stats(sort_columns).print_stats(amount)

    return print_stats
