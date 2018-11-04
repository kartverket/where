"""Plot residuals to console

Description:
------------

Example:
--------

"""

# External library imports
import hipsterplot

# Midgard imports
from midgard.dev import console

# Where imports
from where.lib import log
from where.lib import plugins


@plugins.register
def plot_residuals(dset):
    """Plot residuals

    Args:
        dset:   Dataset, information about model run.
    """
    log.out(f"Residuals at stage {dset.vars['stage']}")
    hipsterplot.plot(
        x_vals=dset.time.utc.datetime, y_vals=dset.residual, num_x_chars=console.columns() - 12, num_y_chars=20
    )
