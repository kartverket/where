"""A Dataset position field

"""

# Midgard imports
from midgard.data.fieldtypes.position import PositionField as MgPositionField
from midgard.dev import plugins

# Where imports
from where.data.position import Position


@plugins.register
class PositionField(MgPositionField):

    _factory = staticmethod(Position)
