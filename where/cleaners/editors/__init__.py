"""Framework for editing data

Description:
------------

Each editor should be defined in a separate .py-file. The function inside the .py-file that should be called
needs to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midgard.dev import plugins

    @plugins.register
    def repair_cycle_slips(dset):
        ...

"""
# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log


def apply_editors(config_key, dset):
    """Apply editors for a given session

    Args:
        config_key (String):  The configuration key listing which editors to apply.
        dset (Dataset):       Dataset containing analysis data.
    """
    prefix = dset.vars["pipeline"]
    editors = config.tech[config_key].list
    log.info(f"Applying editors")
    return plugins.call_all(package_name=__name__, plugins=editors, prefix=prefix, dset=dset)
