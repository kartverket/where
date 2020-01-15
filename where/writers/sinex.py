"""Write normal equations on SINEX format for submission to IVS

Description:
------------

Write normal equations to the SINEX format version 2.02 described in [1]_.

@todo: get solution code ccccc from IVS

References:
-----------

.. [1] SINEX Format https://www.iers.org/IERS/EN/Organization/AnalysisCoordinator/SinexFormat/sinex.html

"""

# Midgard imports
from midgard.dev import plugins
from midgard.files import dependencies

# Where imports
from where.writers import _sinex_blocks_202 as sinex_blocks
from where.lib import config
from where.lib import log

WRITER = __name__.split(".")[-1]


@plugins.register
def write_sinex(dset):
    """Write normal equations of session solution in SINEX format.

    Args:
        dset:  Dataset, data for a model run.
    """
    # Add dependency to sinex_blocks-module
    dependencies.add(sinex_blocks.__file__)

    if config.tech.analysis_status.status.str == "bad":
        log.info("Bad session. Not producing SINEX.")
        return
    with config.files.open("output_sinex", file_vars=dset.vars, mode="wt") as fid:
        sinex = sinex_blocks.SinexBlocks(dset, fid)
        sinex.header_line()
        for block in config.tech[WRITER].blocks.list:
            block_name, *args = block.split(":")
            sinex.write_block(block_name, *args)
        sinex.end_line()
