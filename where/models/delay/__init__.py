"""Framework for calculating delay models

Description:
------------

Each delay model should be defined in a separate .py-file. The function inside the .py-file that should be called need
to be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def vlbi_vacuum_delay(dset):
        ...

The decorated function will be called with a single parameter, ``dset`` which contains a
:class:`~where.data.dataset.Dataset` with data that can be used when calculating the delay.




"""

# Where imports
from where.lib import config
from where.lib import plugins


def calculate(config_key, dset):
    prefix = config.analysis.get("analysis", default="").str
    return plugins.call_all(package_name=__name__, config_key=config_key, prefix=prefix, dset=dset)


def doc(config_key, long_doc=True):
    prefix = config.analysis.get("analysis", default="").str
    return plugins.doc_all(package_name=__name__, config_key=config_key, prefix=prefix, long_doc=long_doc)
