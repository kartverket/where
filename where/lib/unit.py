"""Where library module for handling of SI-unit conversions

Description:
------------

See midgard.math.unit for full documentation

Note that `pint` has a system for defining new units and constants if necessary,
`http://pint.readthedocs.io/en/latest/defining.html`. To use this system, add units to the `units.conf` file in the
`config`-directory.
"""

# Midgard imports
from midgard.math.unit import Unit

# Where imports
from where.lib import files

# Read extra units defined specially for Where
with files.open("units") as fid:
    Unit._ureg.load_definitions(fid)
