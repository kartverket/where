"""Framework for running estimators

Description:
------------

Each estimator should be defined in a separate .py-file. The function inside the .py-file that
should be called need to be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def estimate_least_square(dset):
        ...

The decorated function will be called with a single parameter, ``dset`` which contains a
:class:`~where.data.dataset.Dataset` with data that can be used in the estimation. If the estimator needs the partial
derivatives it should obtain them by itself calling the :func:`where.estimation.partial_vectors` function::

    from where import estimation
    partial_vectors = estimation.partial_vectors(dset, 'estimate', 'estimate_stochastic')



$Revision: 15116 $
$Date: 2018-05-15 17:11:13 +0200 (Tue, 15 May 2018) $
$LastChangedBy: dahmic $

"""
import numpy as np

# Where imports
from where.lib import config
from where.lib import plugins


def get(dset, param_names):
    """Call an ..

    Args:
        dset (Dataset):          Model run data.
        param_names:             Names of parameters to estimate
    """
    constraints = plugins.call_all(
        package_name=__name__, config_key="estimate_constraint", prefix="todo", dset=dset, param_names=param_names
    )

    h = np.concatenate([c[0] for c in constraints.values()])
    sigma = np.concatenate([c[1] for c in constraints.values()])
    import IPython

    IPython.embed()
    return h, sigma
