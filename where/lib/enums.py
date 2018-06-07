"""Definition of Where-specific enumerations

Description:
------------

Custom enumerations used by Where for structured names.


"""

# Standard library imports
import enum

# Make Midgard-enums functions available
from midgard.collections.enums import get_enum, get_value, register_enum  # noqa


#
# ENUMS
#
@register_enum("write_level")
class WriteLevel(enum.IntEnum):
    """Levels used when deciding which fields of a dataset and other information to write to disk"""

    detail = enum.auto()
    analysis = enum.auto()
    operational = enum.auto()


@register_enum("reference_ellipsoid")
class ReferenceEllipsoid(enum.IntEnum):
    """IDs for reference ellipsoids as used by the Sofa library"""

    wgs84 = 1
    grs80 = 2
    wgs72 = 3
