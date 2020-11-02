"""Framework for post-processing data

Description:
------------

Each postprocessor should be defined in a separate .py-file. The function inside the .py-file that should be called
needs to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midgard.dev import plugins

    @plugins.register
    def gnss_linear_combination(dset):
        ...

"""
# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log


def apply_postprocessors(config_key, dset):
    """Apply postprocessors for a given session

    Args:
        config_key (String):  The configuration key listing which postprocessors to apply.
        dset (Dataset):       Dataset containing analysis data.
    """
    prefix = dset.vars["pipeline"]
    postprocessors = config.tech[config_key].list
    log.info(f"Applying postprocessors")
    return plugins.call_all(package_name=__name__, plugins=postprocessors, prefix=prefix, dset=dset)
