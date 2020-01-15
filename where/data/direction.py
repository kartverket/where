"""Array with positions
"""

from typing import Any

# Third party imports
import numpy as np

# Where imports
from where.data import _direction
from where.data._direction import DirectionArray
from where.lib import transformation as trans


def Direction(val=None, ra=None, dec=None, system=None, **dir_args) -> "DirectionArray":
    """Factory for creating PositionArrays for different systems

    See each direction class for exact optional parameters.
    A direction class will always be created in the gcrs system.

    Args:
        val:       Array of right ascension (first column) and declination  (second column) values.
        pos_args:  Additional arguments used to create the DirectionArray.

    Returns:
        Array with positions in the given system.
    """
    if val is None and ra is not None and dec is not None:
        unit_vector = np.array((np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec))).T
        system = "gcrs"
    elif val is not None and system is not None and ra is None and dec is None:
        unit_vector = val
    else:
        raise ValueError(
            "Direction object must be instansiated either with right ascension (ra) and declination (dec)"
            " or a unit direction vector (val) and system"
        )
    return DirectionArray.create(unit_vector, system=system, **dir_args)


#
# Attributes
#
_direction.register_attribute("time")

#
# Position systems
#
@_direction.register_system(convert_to=dict(gcrs=trans.t2g_pos))
class TrsDirection(DirectionArray):

    system = "trs"
    column_names = ("x", "y", "z")
    _units = ("unitless", "unitless", "unitless")


@_direction.register_system(convert_to=dict(trs=trans.g2t_pos))
class GcrsDirection(DirectionArray):

    system = "gcrs"
    column_names = ("x", "y", "z")
    _units = ("unitless", "unitless", "unitless")
