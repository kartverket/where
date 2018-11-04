"""Set up a plug-in architecture for Where

Description:
------------

In order to be able to add models, parsers, data sources etc without needing to hardcode names, but rather pick them
from configuration files, we use a simple plug-in architecture. The plug-in mechanism is based on the different
plug-ins registering themselves using the :func:`register` decorator::

    from where.lib import plugins

    @plugins.register
    def simple_model(rundate, tech, dset):
        ...

Plug-ins are registered based on the name of the module (file) they are defined in, as well as the package (directory)
which contains them. Typically all plug-ins of a given type are collected in a package, e.g. models, techniques,
parsers, etc. To list all plug-ins in a package use :func:`list_all`::

    > from where.lib import plugins
    > plugins.list_all('where.models')
    ['model_one', 'model_three', 'model_two']

If the optional parameter `config_key` is given, then only plug-ins listed in the corresponding section in the current
configuration file is listed. For instance, if the configuration file contains a line saying

    ham_models = model_three, model_one

then we can list only the `ham_models` as follows::

    > from where.lib import plugins
    > plugins.list_all('where.models', config_key='ham_models')
    ['model_one', 'model_three']

Note that the plug-ins by default are sorted alphabetically.

To run the plug-ins, use either :func:`call_all` or :func:`call_one`. The former calls all plug-ins and returns a
dictionary containing the result from each plug-in. As with :func:`list_all` the optional parameter `config_key` may be
given::

    > from where.lib import plugins
    > plugins.call_all('where.models', config_key='ham_models', arg_to_plugin='hello')
    {'model_three': <result from model_three>, 'model_one': <result from model_one>}

Arguments to the plug-ins should be passed as named arguments to :func:`call_all`.

Similarly, one plug-in may be called explicitly using :func:`call_one`::

    > from where.lib import plugins
    > plugins.call_one('where.models', plugin_name='model_one', arg_to_plugin='hello')
    <result from model_one>

There may be more than one function in each plug-in that is decorated by :func:`register`. In this case, the default
behaviour is that only the first function will be called. To call the other registered functions one should use the
:func:`list_parts` function to get a list of these functions and call them explicitly using the `part` optional
parameter to :func:`call_one`::

    > from where.lib import plugins
    > plugins.list_parts('where.techniques', plugin_name='vlbi')
    ['read', 'edit', 'calculate', 'estimate', 'write_result'])
    > for part in plugins.list_parts('where.techniques', plugin_name='vlbi'):
    ...   plugins.call_one('where.techniques', plugin_name='vlbi', part=part, ...)

"""
# Standard library imports
from collections import namedtuple
import importlib
import pathlib
import re
import sys

# Where imports
from where.lib import cache
from where.lib import config
from where.lib import dependencies
from where.lib import exceptions
from where.lib import log
from where.lib.timer import timer


# The _PLUGINS-dict is populated by the :func:`register` decorator in each module.
_PLUGINS = dict(__loaded__=dict())

# Simple structure containing information about a plug-in
Plugin = namedtuple("Plugin", ["name", "function", "file_path", "sort_value"])
Plugin.__doc__ = """Information about a plug-in

    Args:
        name (String):        Name of the plug-in.
        function (Function):  The plug-in.
        file_path (Path):     Path to the source code of the plug-in, used to add the source as a dependency.
        sort_value (Number):  Value used when sorting plug-ins in order to control the order they are called.
    """


#
# REGISTER PLUG-INS
#
def register(func, name=None, sort_value=0):
    """Register a plug-in

    Plug-ins are registered based on the name of the module (file) they are defined in, as well as the package
    (directory) which contains them. Typically all plug-ins of a given type are collected in a package, e.g. models,
    techniques, parsers, etc. The path to the source code file is also stored. This is used to be able to add the
    source code as a dependency file when the plug-in is called.

    If `name` is given, the plug-in is registered based on this name instead of the name of the module. The name of the
    module is still registered as a part that can be used to distinguish between similar plug-ins in different files
    (see for instance how `session` is used in `where.techniques`).

    Args:
        func (Function):       The function that is being registered.
        name (String):         Alternative name of plug-in. Used by `register_named`.
        sort_value (Number):   The value used when sorting plug-ins. Used by `register_ordered`.

    Returns:
        Function: The function that is being registered.
    """
    # Get information from the function being registered
    package_name, _, plugin_name = func.__module__.rpartition(".")
    file_path = pathlib.Path(sys.modules[func.__module__].__file__)

    # Store Plugin-object in _PLUGINS dictionary
    plugin_info = _PLUGINS.setdefault(package_name, dict()).setdefault(plugin_name, dict())
    if name is None:
        name = func.__name__  # Name of function is used as default name
        plugin_info.setdefault("__parts__", list()).append(name)  # Only unnamed parts are added to list

    plugin = Plugin("{}.{}".format(plugin_name, name), func, file_path, sort_value)
    plugin_info[name] = plugin
    log.debug("Registering {} as a {}-plugin from {}", plugin.name, package_name, plugin.file_path)

    # Add first registered unnamed part as default
    if "__parts__" in plugin_info:
        plugin_info["__default__"] = plugin_info[plugin_info["__parts__"][0]]

    return func


def register_named(name):
    """Register a named plug-in

    This allows for overriding the name used to register the plug-in. See `register` for more details.

    Args:
        name (String):   Name used for plug-in instead of module name.

    Returns:
        Decorator: Decorator that registers a function.
    """

    def register_decorator(func):
        return register(func, name=name)

    return register_decorator


def register_ordered(sort_value):
    """Register a plug-in with a specific sort order

    The sort value should be a number. Lower numbers are sorted first, higher numbers last. Plug-ins without an
    explicit sort_order gets the sort value of 0.

    Args:
        sort_value (Number):   The value used when sorting plug-ins.

    Returns:
        Decorator: Decorator that registers a function.
    """

    def register_decorator(func):
        return register(func, sort_value=sort_value)

    return register_decorator


#
# CALL PLUG-INS
#
def call_one(
    package_name, plugin_name, part=None, prefix=None, logger=log.time, use_timer=True, do_report=True, **kwargs
):
    """Call one plug-in

    If the plug-in is not part of the package an UnknownPluginError is raised.

    If there are several functions registered in a plug-in and `part` is not specified, then the first function
    registered in the plug-in will be called.

    The file containing the source code of the plug-in is added to the list of dependencies.

    Args:
        package_name (String):  Name of package containing plug-ins.
        plugin_name (String):   Name of the plug-in, i.e. the module containing the plug-in.
        part (String):          Name of function to call within the plug-in (optional).
        prefix (String):        Prefix of the plug-in name, used if the plug-in name is unknown (optional).
        logger (Function):      Logger from the lib.log package specifying the level of logging to be used (optional).
        use_timer (Boolean):    Whether to time and log the call to the plug-in (optional).
        do_report (Boolean):    Whether to add the call to the plug-in to the report (optional).
        kwargs:                 Named arguments passed on to the plug-in.

    Returns:
        Return value of the plug-in.
    """
    # Get Plugin-object
    plugin_name = load_one(package_name, plugin_name, prefix=prefix)
    part = "__default__" if part is None else part
    try:
        plugin = _PLUGINS[package_name][plugin_name][part]
    except KeyError:
        raise exceptions.UnknownPluginError(
            "Plugin '{}' not found for '{}' in '{}'" "".format(part, plugin_name, package_name)
        ) from None

    # Add plug-in to report
    if do_report:
        from where.reports import report

        code_kwargs = kwargs.copy()
        if "dset" in code_kwargs:
            code_kwargs["dset"] = code_kwargs["dset"].repr
        report.add(
            package_name,
            __plugin__=plugin.name,
            __doc__=plugin.function.__doc__,
            __text__="TODO",
            __code__="kwargs = {}\n{} = plugins.call_one('{}', '{}', part='{}', **kwargs)"
            "".format(code_kwargs, plugin_name, package_name, plugin_name, part),
            **kwargs,
        )

    # Call plug-in
    dependencies.add(plugin.file_path)
    if logger:
        logger("Start {} in {}", plugin.name, package_name)
        time_logger = log.lowest(logger, log.time) if use_timer else None
    else:
        time_logger = None
    with timer("Finish {} ({}) in".format(plugin.name, package_name), logger=time_logger):
        return plugin.function(**kwargs)


def call_all(
    package_name, plugins=None, config_key=None, prefix=None, logger=log.time, use_timer=True, do_report=True, **kwargs
):
    """Call all plug-ins in a package

    Either `plugins` or `config_key` may be specified. If `plugins` is given, it should be a list of names of plug-ins.
    If `config_key` is given, only plug-ins listed in the config file with the given key are called. If a plug-in
    listed in the `plugins`-list or in the config file does not exist, an UnknownPluginError is raised.

    If neither `plugins` nor `config_key` is given, all available plugins will be called. Do note, however, that this
    will import all python files in the package.

    Args:
        package_name (String):  Name of package containing plug-ins.
        plugins (Tuple):        List of plug-ins that should be used (optional).
        config_key (String):    Key in config file containing list of plug-ins that should be used (optional).
        prefix (String):        Prefix of the plug-in names, used if any of the plug-in names are unknown (optional).
        logger (Function):      Logger from the lib.log package specifying the level of logging to be used (optional).
        use_timer (Boolean):    Whether to time and log the calls to the plug-ins (optional).
        do_report (Boolean):    Whether to add calls to the plug-ins to the report (optional).
        kwargs:                 Named arguments passed on to all the plug-ins.

    Returns:
        Dict: Dictionary of all results from the plug-ins.
    """
    plugin_names = list_all(package_name, plugins=plugins, config_key=config_key, prefix=prefix)
    return {
        p: call_one(package_name, p, logger=logger, use_timer=use_timer, do_report=do_report, **kwargs)
        for p in plugin_names
    }


#
# GET DOCUMENTATION FOR PLUG-INS
#
def doc_one(package_name, plugin_name, part=None, prefix=None, long_doc=True, include_details=False):
    """Document one plug-in

    If the plug-in is not part of the package an UnknownPluginError is raised.

    If there are several functions registered in a plug-in and `part` is not specified, then the first function
    registered in the plug-in will be documented.

    Args:
        package_name (String):     Name of package containing plug-ins.
        plugin_name (String):      Name of the plug-in, i.e. the module containing the plug-in.
        part (String):             Name of function to call within the plug-in (optional).
        prefix (String):           Prefix of the plug-in name, used if the plug-in name is unknown (optional).
        long_doc (Boolean):        Whether to return the long doc-string or the short one-line string (optional).
        include_details (Boolean): Whether to include development details like parameters and return values (optional).

    Returns:
        String: Documentation of the plug-in.
    """
    # Get Plugin-object and pick out doc-string
    plugin_name = load_one(package_name, plugin_name, prefix=prefix)
    part = "__default__" if part is None else part
    try:
        plugin = _PLUGINS[package_name][plugin_name][part]
    except KeyError:
        raise exceptions.UnknownPluginError(
            "Plugin '{}' not found for '{}' in '{}'" "".format(part, plugin_name, package_name)
        ) from None
    doc = plugin.function.__doc__ if plugin.function.__doc__ else ""

    if long_doc:
        # Strip short description and indentation
        lines = [d.strip() for d in "\n\n".join(doc.split("\n\n")[1:]).split("\n")]

        # Stop before Args:, Returns: etc if details should not be included
        idx_args = len(lines)
        if not include_details:
            re_args = re.compile("(Args:|Returns:|Details:|Examples?:|Attributes:)$")
            try:
                idx_args = [re_args.match(l) is not None for l in lines].index(True)
            except ValueError:
                pass
        return "\n".join(lines[:idx_args]).strip()

    else:
        # Return short description
        return doc.split("\n\n")[0].replace("\n", " ").strip()


def doc_all(package_name, plugins=None, config_key=None, prefix=None, long_doc=True, include_details=False):
    """Call all plug-ins in a package

    Either `plugins` or `config_key` may be specified. If `plugins` is given, it should be a list of names of plug-ins.
    If `config_key` is given, only plug-ins listed in the config file with the given key are called. If a plug-in
    listed in the `plugins`-list or in the config file does not exist, an UnknownPluginError is raised.

    If neither `plugins` nor `config_key` is given, all available plugins will be called. Do note, however, that this
    will import all python files in the package.

    Args:
        package_name (String):     Name of package containing plug-ins.
        plugins (Tuple):           List of plug-ins that should be used (optional).
        config_key (String):       Key in config file containing list of plug-ins that should be used (optional).
        prefix (String):           Prefix of the plug-in names, used if any of the plug-ins are unknown (optional).
        long_doc (Boolean):        Whether to return the long doc-string or the short one-line string (optional).
        include_details (Boolean): Whether to include development details like parameters and return values (optional).

    Returns:
        Dict: Dictionary of all results from the plug-ins.
    """
    plugin_names = list_all(package_name, plugins=plugins, config_key=config_key, prefix=prefix)
    return {p: doc_one(package_name, p, long_doc=long_doc, include_details=include_details) for p in plugin_names}


#
# LIST AVAILABLE PLUG-INS
#
@cache.function
def list_all(package_name, plugins=None, config_key=None, prefix=None):
    """Load plug-ins in a package

    Either `plugins` or `config_key` may be specified. If `plugins` is given, it should be a list of names of plug-ins.
    If `config_key` is given, only plug-ins listed in the config file with the given key are returned. If a plug-in
    listed in the `plugins`-list or in the config file does not exist, an UnknownPluginError is raised.

    If neither `plugins` nor `config_key` is given, all available plugins will be listed. Do note, however, that this
    will import all python files in the package.

    Args:
        package_name (String):  Name of package containing plug-ins.
        plugins (Tuple):        List of plug-ins that should be used (optional).
        config_key (String):    Key in config file containing list of plug-ins that should be used (optional).
        prefix (String):        Prefix of the plug-in names, used if any of the plug-in names are unknown (optional).

    Returns:
        List: List of strings with names of plug-ins.
    """
    if plugins is not None and config_key is not None:
        raise exceptions.InitializationError("Either 'plugins' or 'config_key' should be specified, not both")

    # Figure out names of plug-ins
    if plugins is not None:
        plugin_names = plugins
    elif config_key is not None:
        plugin_names = config.tech[config_key].tuple
    else:
        _import_all(package_name)
        plugin_names = _PLUGINS[package_name].keys()

    # Load each plug-in and return them in sort order
    def _sort_value(plugin):
        """Pick out sort_value of plugin"""
        return _PLUGINS[package_name][plugin]["__default__"].sort_value, plugin

    return sorted((load_one(package_name, p, prefix=prefix) for p in plugin_names), key=_sort_value)


@cache.function
def list_parts(package_name, plugin_name, prefix=None):
    """List all parts of one plug-in

    Args:
        package_name (String):  Name of package containing plug-ins.
        plugin_name (String):   Name of the plug-in.
        prefix (String):        Prefix of the plug-in name, used if the plug-in name is unknown (optional).

    Returns:
        List: Strings with names of parts.
    """
    plugin_name = load_one(package_name, plugin_name, prefix=prefix)
    return _PLUGINS[package_name][plugin_name].get("__parts__", list())


#
# LOAD PLUG-INS
#
@cache.function
def load_one(package_name, plugin_name, prefix=None):
    """Load one plug-in from a package

    First tries to load the plugin with the given name. If that fails, it tries to load {prefix}_{plugin_name} instead.

    Args:
        package_name (String):  Name of package containing plug-ins.
        plugin_name (String):   Name of the plug-in (module).
        prefix (String):        Prefix of the plug-in name, used if the plug-in name is unknown (optional).

    Returns:
        String: Actual name of plug-in (with or without prefix).
    """
    if plugin_name not in _PLUGINS.get(package_name, dict()):
        try:
            _import_one(package_name, plugin_name)
        except exceptions.UnknownPluginError as import_err:
            if prefix:
                try:
                    return load_one(package_name, f"{prefix}_{plugin_name}")
                except exceptions.UnknownPluginError:
                    pass

            raise import_err from None

    return plugin_name


def _import_one(package_name, plugin_name):
    """Import a plugin from a package

    This is essentially just a regular python import. As the module is imported, the _PLUGINS-dict will be populated by
    @register decorated functions in the file.

    Args:
        package_name (String):  Name of package containing plug-ins.
        plugin_name (String):   Name of the plug-in (module).
    """
    try:
        importlib.import_module(package_name + "." + plugin_name)
    except ImportError:
        raise exceptions.UnknownPluginError(
            "Plug-in '{}' not found in package '{}'" "".format(plugin_name, package_name)
        ) from None


def _import_all(package_name):
    """Import the relevant .py-files in the given package directory

    As each file is imported, the _PLUGINS-dict will be populated by @register decorated functions in the files.

    Args:
        package_name (String):  Name of package containing plug-ins.
    """
    # Figure out the directory of the package by importing it
    try:
        package = importlib.import_module(package_name)
    except ImportError:
        raise exceptions.UnknownPluginError("Plug-in package '{}' not found".format(package_name)) from None

    # Import all .py files in the given directory
    directory = pathlib.Path(package.__file__).parent
    for file_path in directory.glob("*.py"):
        plugin_name = file_path.stem
        if not plugin_name.startswith("_"):
            _import_one(package_name, plugin_name)
