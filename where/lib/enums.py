"""Definition of Where-specific enumerations

Description:
------------

Custom enumerations used by Where for structured names.


"""

# Standard library imports
import colorama
import enum

# Make Midgard-enums functions available
from midgard.collections.enums import get_enum, get_value, register_enum  # noqa


#
# ENUMS
#
@register_enum("reference_ellipsoid")
class ReferenceEllipsoid(enum.IntEnum):
    """IDs for reference ellipsoids as used by the Sofa library"""

    wgs84 = 1
    grs80 = 2
    wgs72 = 3


@register_enum("log_level")
class LogLevel(int, enum.Enum):
    """Levels used when deciding how much log output to show"""

    all = enum.auto()
    debug = enum.auto()
    time = enum.auto()
    dev = enum.auto()
    info = enum.auto()
    out = enum.auto()
    warn = enum.auto()
    check = enum.auto()
    error = enum.auto()
    fatal = enum.auto()
    none = enum.auto()


@register_enum("log_color")
class LogColor(str, enum.Enum):
    """Colors used when logging"""

    dev = (colorama.Fore.BLUE,)
    time = (colorama.Fore.WHITE,)
    out = (colorama.Style.BRIGHT,)
    check = (colorama.Style.BRIGHT + colorama.Fore.YELLOW,)
    warn = colorama.Fore.YELLOW
    error = colorama.Fore.RED
    fatal = colorama.Style.BRIGHT + colorama.Fore.RED
