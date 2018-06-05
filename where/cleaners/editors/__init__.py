"""Framework for editing data

Description:
------------

Each editor should be defined in a separate .py-file. The function inside the .py-file that should be called
needs to be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def repair_cycle_slips(dset):
        ...

"""

# Where imports
from where.lib import config
from where.lib import log
from where.lib import plugins


def apply_editors(config_key, dset):
    """Apply editors for a given session

    Args:
        config_key (String):  The configuration key listing which editors to apply.
        dset (Dataset):       Dataset containing analysis data.
    """
    prefix = config.analysis.get("analysis", default="").str
    editors = plugins.list_all(package_name=__name__, config_key=config_key, prefix=prefix)
    log.info(f"Applying editors {', '.join(editors)}")
    return plugins.call_all(package_name=__name__, config_key=config_key, prefix=prefix, dset=dset)
