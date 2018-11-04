"""Where library module defining an assortment of constants

Example:

    >>> from where.lib import constant
    >>> print('The speed of light is {:0.2f}'.format(constant.c))
    The speed of light is 299792458.00


Description:
------------

This module provides constants that are used within the Where project. The actual constants are defined in the
`constants.conf` file (see the file list for location). See that file for references and for adding or changing
constants.

The constants are stored as module variables so they can be used simply as `constant.c` as in the example above. Some
models use particular values for constants that are different from the conventional ones. This is handled by the source
parameter. For instance, the EGM 2008 gravity field is calculated with a value for GM different from the IERS
Conventions value, using::

    constant.get('GM', source='egm_2008')

instead of simply `constant.GM`.


Todo:
-----

Rewrite as a class instead of a module, to have somewhat cleaner code (and be more consistent with things like
lib.unit).

"""
# Standard library imports
from configparser import ConfigParser

# Where imports
from where.lib import files
from where.lib import log

# Where constants as a ConfigParser, loaded from file.
# TODO: Use midgard.config.Configuration
_CONSTANTS = ConfigParser()


def get(constant, source="default"):
    """Get the value of one constant

    Note that if you need the default value of the constant (from the default source) it is typically better to read it
    as a property. That is, `constant.c` is preferred to `constant.get('c')`.

    Args:
        constant:   Name of the constant.
        source:     Source from which the constant is defined.

    Returns:
        float: Value of constant.
    """
    try:
        return float(_CONSTANTS[constant][source])
    except KeyError:
        raise KeyError("Constant '{}' is not defined by source '{}'".format(constant, source)) from None


def use_source(source="default"):
    """Use a given source for constants

    This overwrites all values of constants with values from the given source. That is, the values given when reading
    constants as properties (e.g. `constant.c`) will be taken from source.

    Take care when using a different source than default that the source is reset to default. You might want to use
    code similar to the following::

        constant.use_source('egm_2008')
        do_something_awesome()
        constant.use_source('default')

    Args:
        source:   Source from which the constant is defined.
    """
    for constant in _CONSTANTS.sections():
        value = _CONSTANTS[constant].get(source)
        if value is not None:
            globals()[constant] = float(value)


def _read_constants_from_file():
    """Read constants from the constants file

    The constants are read from the constants file defined in the filelist, and stored in a ConfigParser for use in the
    :func:`get` and :func:`use_source` functions later.
    """
    constant_path = files.path("constants")
    log.debug("Use Where constants from {}", constant_path)
    _CONSTANTS.read(constant_path)

    # Set up aliases
    aliased_constants = [c for c in _CONSTANTS.sections() if _CONSTANTS.has_option(c, "_aliases")]
    for constant in aliased_constants:
        aliases = _CONSTANTS[constant].get("_aliases").replace(",", " ").split()
        for alias in aliases:
            if alias in _CONSTANTS:
                raise KeyError(
                    "Constant '{}' is already defined. The name '{}' can not be added as an alias for '{}'"
                    " in {}".format(alias, alias, constant, constant_path)
                )
            _CONSTANTS[alias] = _CONSTANTS[constant]


def load():
    """Read the constants and store as global variables in the module

    The default values for the constants are added as global members of the module, and can be accessed for instance as
    `constant.c`.
    """
    _read_constants_from_file()
    use_source()


# Read the constants immediately when the module is imported
load()
