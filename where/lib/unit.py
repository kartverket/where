"""Where library module for handling of SI-unit conversions

Description:
------------

This module provides unit conversion constants and functions. The heavy lifting is done by the `pint` package. The
basic usage is as follows::

    from where.lib.unit import unit
    seconds_in_two_weeks = 2 * unit.week2secs    # -> 1209600.0

In general `unit.spam2ham` will give the multiplicative conversion scale between the units `spam` and `ham`. Through
the `pint` package we support a lot of units. See :func:`unit.list()` or
`https://github.com/hgrecco/pint/blob/master/pint/default_en.txt`. Another notation is also available, and might be
necessary for some more complicated conversions::

    seconds_in_two_weeks = 2 * unit('week', 'seconds')
    miles_per_hour_in_meters_per_second = unit('mph', 'meters / sec')

Do note that we support most normal aliases as well as singular and plural forms of the units. For instance can
`second` be represented as `s`, `sec`, `secs` and `seconds`. Prefixes are also handled::

    nanoseconds_in_an_hour = unit.hour2nanosecs
    inches_in_a_kilometer = unit.km2inches

For more complicated conversions (for instance from Celsius to Fahrenheit) one can create custom conversion functions
using :func:`convert`::

    c2f = unit.convert('celsius', 'fahrenheit')
    absolute_zero_in_fahrenheit = c2f(-273.15)

For convenience, this can also be written using the attribute notation as `unit.spam_to_ham(spam_value)`. Then the
previous example simply becomes::

    absolute_zero_in_fahrenheit = unit.celsius_to_fahrenheit(-273.15)

(or even easier `unit.kelvin_to_fahrenheit(0)`).

Finally, we can access the unit/quantity system of `pint` by using the name of a unit by itself, e.g. `unit.spam`. For
instance::

    distance = 42 * unit.km
    time = 31 * unit('minutes')
    speed = distance / time
    speed.to(unit.mph)          # -> <Quantity(50.511464659292955, 'mph')>
    speed.to_base_units()       # -> <Quantity(22.580645161290324, 'meter / second')>

However, using the full unit system adds some overhead so we should be careful in using it in heavy calculations.

Note that `pint` has a system for defining new units and constants if necessary,
`http://pint.readthedocs.io/en/latest/defining.html`. To use this system, add units to the `units.conf` file in the
`config`-directory.





"""

# External library imports
import numpy as np
import pint

# Where imports
from where.lib import cache
from where.lib import files
from where.lib.exceptions import UnitError

# The _UNITS-dict is used to keep track of units values returned by functions and methods
_UNITS = dict()


class _convert_units(type):
    """A meta-class that does the parsing of units

    The meta-class is used for convenience. It allows us to use the `unit`-class without instantiating it. That is, we
    can write `unit.km2m` instead of `unit().km2m`.
    """

    ureg = pint.UnitRegistry()

    with files.open("units") as fid:
        ureg.load_definitions(fid)

    @cache.function
    def __call__(cls, from_unit, to_unit=None):
        """Calculate the conversion scale between from_unit and to_unit

        If `to_unit` is not given, then `from_unit` is interpreted as a constant which is converted to base units
        (meters, seconds, etc) and returned.

        Args:
            from_unit (String):   The unit to convert from
            to_unit (String):     The unit to convert to

        Returns:
            Float:  Scale to multiply by to convert from from_unit to to_unit
        """
        if to_unit is None:
            return cls.ureg(from_unit)
        else:
            return cls.ureg(from_unit).to(to_unit).magnitude

    def __getattr__(cls, key):
        """Simplify notation for converting between units

        This makes it possible to type `unit.km2m` instead of `unit('km', 'm')`. We split on the character `2`
        (pronounced "to"), and pass the result on to :func:`__call__` to do the conversion. If a `2` is not found, we
        check if we can split on '_to_' instead, if so it is interpreted as a conversion function and is handed of to
        :func:`convert`. Finally, if no split is done, the attribute is interpreted as a simple unit.

        Note that if you need a unit whose name contains a '2' (or '_to_') you need to use the notation
        `unit('foot_H2O', 'pascal'). Similarly, more complex units need the same notation, e.g. `unit('meters per
        second ** 2')`.

        Args:
            key (String):   The key (name) of the attribute to the class. Interpreted as units

        Returns:
            Float:  Scale to multiply by to perform the unit conversion
        """
        if "2" in key:
            from_unit, to_unit = key.split("2", maxsplit=1)
            return cls(from_unit, to_unit)
        elif "_to_" in key:
            from_unit, to_unit = key.split("_to_", maxsplit=1)
            return cls.function(from_unit, to_unit)
        else:
            return cls(key)

    def function(cls, from_unit, to_unit):
        """Create a conversion function

        This is necessary for unit conversions that are not simple multiplications. The usual example is temperature
        conversions for instance from Celsius to Fahrenheit.

        Args:
            from_unit (String):   The unit to convert from
            to_unit (String):     The unit to convert to

        Returns:
            Function:  Conversion function that converts from from_unit to to_unit
        """
        return lambda value: cls.ureg.Quantity(value, cls.ureg(from_unit)).to(cls.ureg(to_unit)).magnitude

    def register(cls, unit):
        """Register unit of a function/method/property

        This method should be used as a decorator on the function/method/property, and specify the unit of the value
        returned by that function/method/property. For instance

            @property
            @unit.register('meter')
            def calculate_delay(...):
                return delay_in_meters

        Units registered with this decorator can be used by the functions returned by the `unit_func_factory`,
        `convert_func_factory` and `factor_func_factory`.

        Args:
            unit (String):  Name of unit.

        Returns:
            Function:  Decorator that registers the unit.
        """

        def register_decorator(func):
            """Register unit of func in _UNITS-dictionary"""
            module_name = func.__module__
            func_name = func.__name__
            _UNITS.setdefault(module_name, dict())[func_name] = unit

            return func

        return register_decorator

    @staticmethod
    def _get_unit(module_name, func_name):
        """Get registered unit of function/method/property

        Outside code should use the `unit_factory` to get registered units.

        Args:
            module_name (String):   Name of module containing function/method/property.
            func_name (String):     Name of function/method/property with registered unit.

        Returns:
            String:  Name of unit.
        """
        units = _UNITS.get(module_name, dict())
        try:
            return units[func_name]
        except KeyError:
            raise UnitError("No unit is registered for '{}' in {}".format(func_name, module_name)) from None

    def unit_factory(cls, module_name):
        """Provide a function that can get registered units of functions/methods/properties

        The function checks for units registered with the unit.register-decorator. It can for instance be added to a
        class as follows:

            unit = staticmethod(unit.unit_func_factory(__name__))

        Args:
            module_name (String):   Name of module as returned by `__name__`.

        Returns:
            Function:  Function that gets unit of values returned by functions.
        """

        def unit(func_name):
            """Unit of value returned by function/method/property

            Args:
                func_name (String):  Name of function/method/property.

            Returns:
                String:  Name of unit.
            """
            return cls._get_unit(module_name, func_name)

        return unit

    def convert_factory(cls, module_name):
        """Provide a function that can convert values of properties to a given unit

        The function checks for units registered with the unit.register-decorator. It can for instance be added to a
        class as follows:

            convert_to = unit.convert_property_factory(__name__)

        Note that unlike the other factories, this one only works for properties.

        Args:
            module_name (String):   Name of module as returned by `__name__`.

        Returns:
            Function:  Function that converts values of properties.
        """

        def convert(self, property_name, to_unit):
            """Convert value of property to another unit

            Args:
                property_name (String):  Name of property.
                to_unit (String):        Name of other unit

            Returns:
                Numeric scalar or array:  Values of property converted to other unit.
            """
            from_unit = cls._get_unit(module_name, property_name)
            factor = cls(from_unit, to_unit)
            return getattr(self, property_name) * factor

        return convert

    def factor_factory(cls, module_name):
        """Provide a function that calculates conversion factor to another unit

        The function finds conversion factors for units registered with the unit.register-decorator. It can for
        instance be added to a class as follows:

            unit_factor = staticmethod(unit.factor_factory(__name__))

        Args:
            module_name (String):   Name of module as returned by `__name__`.

        Returns:
            Function:  Function that calculates conversion factor to another unit.
        """

        def factor(func_name, to_unit):
            """Conversion factor between unit of function/method/property and another unit

            Args:
                func_name (String):  Name of function/method/property.
                to_unit (String):    Name of other unit.

            Returns:
                Float:  Conversion factor.
            """
            from_unit = cls._get_unit(module_name, func_name)
            return cls(from_unit, to_unit)

        return factor

    def units_dict(cls, module_name):
        """Dictionary of units registered on a module

        Add a sub-dictionary if the module name is unknown, to set up a reference in case units are registered later.

        Returns:
            Dictionary:  Units registered on a module.
        """
        _UNITS.setdefault(module_name, dict())
        return _UNITS[module_name]

    @property
    def names(cls):
        """List available units and constants

        The list of available units contains aliases (for instance s, sec, second), but not plural forms (secs,
        seconds) or possible prefixes (milliseconds, usec, ms).

        Returns:
            List of strings: Available units and constants
        """
        return dir(cls.ureg)

    #
    # Conversion routines not defined by pint
    #
    def rad_to_dms(cls, radians):
        """Converts radians to degrees, minutes and seconds

        Args:
            radians (Float):  angle(s) in radians

        Returns:
            Tuple of Floats:  degrees, minutes, seconds

        Examples:
            >>> unit.rad_to_dms(1.04570587646256)
            (59.0, 54.0, 52.3200000000179)
            >>> unit.rad_to_dms(-0.2196050301753194)
            (-12.0, 34.0, 56.78900000000468)
            >>> unit.rad_to_dms(-0.005817642339636369)
            (-0.0, 19.0, 59.974869999999925)
        """
        sign = np.sign(radians)
        degrees = abs(radians) * cls.radians2degrees
        minutes = (degrees % 1) * cls.hour2minutes
        seconds = (minutes % 1) * cls.minute2seconds

        return sign * np.floor(degrees), np.floor(minutes), seconds

    def dms_to_rad(cls, degrees, minutes, seconds):
        """Convert degrees, minutes and seconds to radians

        The sign of degrees will be used. In this case, be careful that the sign
        of +0 or -0 is correctly passed on. That is, degrees must be specified as a float, not an
        int.

        Args:
            degrees:   Degrees as float (including sign) or array of floats
            minutes:   Minutes as int/float or array of ints/floats
            seconds:   Seconds as float or array of floats

        Returns:
            Float/Array: Given degrees, minutes and seconds as radians.

        Examples:
            >>> unit.dms_to_rad(59, 54, 52.32)
            1.04570587646256
            >>> unit.dms_to_rad(-12.0, 34, 56.789)
            -0.21960503017531938
            >>> unit.dms_to_rad(-0.0, 19, 59.974870)
            -0.005817642339636369
        """
        sign = np.copysign(1, degrees)
        return (
            sign * (np.abs(degrees) + minutes * cls.minutes2hours + seconds * cls.seconds2hours) * cls.degrees2radians
        )

    def hms_to_rad(cls, hours, minutes, seconds):
        """Convert hours, minutes and seconds to radians

        Args:
            hours:     Hours as int or array of ints
            minutes:   Minutes as int or or array of ints
            seconds:   Seconds as float or or array of floats

        Returns:
            Float: Given hours, minutes and seconds as radians.

        Examples:
            >>> unit.hms_to_rad(17, 7, 17.753427)
            4.482423920139868
            >>> unit.hms_to_rad('12', '0', '0.00')
            3.1415926535897936
            >>> unit.hms_to_rad(-12, 34, 56.789)
            Traceback (most recent call last):
            ValueError: hours must be non-negative
        """
        return 15 * cls.dms_to_rad(hours, minutes, seconds)


class unit(metaclass=_convert_units):
    """Unit converter

    The implementation of the unit conversion is done in the `_convert_units`-metaclass.
    """

    # Make pint exceptions available
    from pint.errors import DimensionalityError
