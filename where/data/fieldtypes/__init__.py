"""Field types that can be used by Dataset

"""
# Midgard imports
from midgard.data.fieldtypes import names, function  # noqa
from midgard.data import fieldtypes as mg_fieldtypes
from midgard.dev import plugins

# Add Where fieldtypes to Midgard fieldtypes
plugins.add_alias(mg_fieldtypes.__name__, __name__)
