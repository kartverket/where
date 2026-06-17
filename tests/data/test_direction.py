# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import exceptions 

# Where imports
from where.data import direction

def test_pos_conversions(src_dir):
    systems = direction.Direction.SYSTEMS

    print(f"Testing systems {systems}")
    for system in systems:
        try:
            converted_src_dir = getattr(getattr(src_dir, system), src_dir.system)
            assert np.allclose(np.asarray(src_dir), np.asarray(converted_src_dir), equal_nan=True)
            print(f"src_dir.{system} == src_dir.{system}.{src_dir.system} OK")
        except exceptions.UnknownConversionError:
            print(f"Conversion from {src_dir.system} to {system} is not defined")