"""Framework for calculating models

Description:
------------

Each model should be defined in a separate .py-file inside one of the model directories, delay, site, orbit and
satellite. The function inside the .py-file that should be called need to be decorated with the
:func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def vlbi_vacuum_delay(dset):
        ...

See the corresponding :file:`__init__.py`-files for more information.


See also:
---------

* :mod:`where.models.delay`
* :mod:`where.models.site`

"""

# Make estimators and parameters functions available
from where.estimation.estimators import *  # noqa
from where.estimation.parameters import *  # noqa
from where.estimation.outliers import *  # noqa
