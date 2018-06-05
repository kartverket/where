"""Write normal equations on SINEX format for submission to IVS

Description:
------------

Write normal equations to the SINEX format version 2.02 described in [1]_.

@todo: get solution code ccccc from IVS

References:
-----------

.. [1] SINEX Format https://www.iers.org/IERS/EN/Organization/AnalysisCoordinator/SinexFormat/sinex.html

"""
# Where imports
from where.writers import sinex_blocks
from where.lib import dependencies
from where.lib import files
from where.lib import plugins
from where.lib import config
from where.lib import log


@plugins.register
def normal_equations(dset):
    """Write normal equations of session solution in SINEX format.

    Args:
        dset:  Dataset, data for a model run.
    """
    # Add dependency to sinex_blocks-module
    dependencies.add(sinex_blocks.__file__)

    if config.tech.analysis_status.status.str == "bad":
        log.info("Bad session. Not producing SINEX.")
        return
    with files.open("output_vlbi_sinex", file_vars=dset.vars, mode="wt") as fid:
        sinex = sinex_blocks.SinexBlocks(dset, fid)
        sinex.header_line()
        sinex.file_reference()
        sinex.file_comment()
        sinex.input_acknowledgements()
        sinex.nutation_data()
        sinex.precession_data()
        sinex.source_id()
        sinex.site_id()
        sinex.site_eccentricity()
        sinex.solution_epochs()
        sinex.solution_statistics()
        sinex.solution_estimate()
        sinex.solution_apriori()
        sinex.solution_normal_equation_vector()
        sinex.solution_normal_equation_matrix("U")
        sinex.end_line()
