"""Framework for cleaning data

Description:
------------

Each data cleaner should be defined in a one of two directories:

+ `editors`  - Editors can add new fields to the dataset.
+ `removers` - These cleaners only remove observations.

"""

# Make the apply-functions in subpackages available
from where.cleaners.editors import apply_editors  # noqa
from where.cleaners.removers import apply_removers  # noqa
from where.cleaners.removers import apply_remover  # noqa

# Do not support * imports
__all__ = []
