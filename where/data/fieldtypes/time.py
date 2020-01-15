"""A Dataset time field

"""

# Midgard imports
from midgard.data.fieldtypes.time import TimeField as MgTimeField
from midgard.dev import plugins

# Where imports
from where.data.time import Time


@plugins.register
class TimeField(MgTimeField):

    _factory = staticmethod(Time)
