"""Framework for cleaning data

Description:
------------

Each data cleaner should be defined in a one of two directories:

+ `editors`  - Editors can add new fields to the dataset.
+ `removers` - These cleaners only remove observations.



$Revision: 14978 $
$Date: 2018-04-30 19:01:11 +0200 (Mon, 30 Apr 2018) $
$LastChangedBy: hjegei $

"""

# Embed configuration functions for easier access
from where.cleaners._config import add_session_config_file, add_session_config, store_session_config  # noqa
from where.cleaners._config import apply_options  # noqa

# Make the apply-functions in subpackages available
from where.cleaners.editors import apply_editors  # noqa
from where.cleaners.removers import apply_removers  # noqa

# Do not support * imports
__all__ = []
