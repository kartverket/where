#!/usr/bin/env python3
"""Read URE control tool CSV output file and convert them to Where dataset

Description:
------------

"""
# Standard library imports
from datetime import date
import subprocess

# External library imports
import numpy as np

# Where imports
from where import data
from where import parsers
from where.lib import log

TECH = "sisre"


def _difference(dset_where, dset_ure_tool):

    if dset_where.num_obs == 0 or dset_ure_tool == 0:
        log.warn(
            f"Nothing to compare. Number of observations are zero at least for one dataset "
            f"(dset_where: {dset_where.num_obs}, dset_ure_tool: {dset_ure_tool.num_obs})."
        )
        return 1

    ddiff = dset_where.difference_with(dset_ure_tool, difference_by=["time", "satellite"])
    ddiff.write_as(tech=TECH, stage="dwhr_ure")


def _generate_where_dset(rundate, nav_type, system, frequencies):
    """Generate Where dataset

    Start Where analysis with Where and read Dataset.
    """
    program = "where"
    session_name = "brdm"

    # Start processing with Where
    where_args = [
        str(rundate.year),
        str(rundate.month),
        str(rundate.day),
        "--" + TECH,
        "--session_name=" + session_name,
        "--brdc_block_nearest_to=transmission_time:positive",
        "--systems=" + system,
        "--frequencies=" + frequencies,
        "--navigation_message_type=" + system + ":" + nav_type,
        "--clock_product=sp3",
        "--weight_factor_elev_mask=5",
        "--profile=grc",
        "-N",
        "-D",
        "-T",
    ]
    log.info("Start {}: {} {}".format(TECH.upper(), program, " ".join(where_args)))
    process = subprocess.run([program] + where_args)
    if process.returncode:
        log.error(f"Where {TECH.upper()} failed with error code {process.returncode} ({' '.join(process.args)})")

    # Write SISRE Dataset
    dset_vars = dict(rundate=rundate, tech="sisre", stage="calculate", session="", dataset_name="", dataset_id="last")
    dset = data.Dataset(**dset_vars)
    dset.add_float(
        "sqrt_a2_c2", val=(np.sqrt(dset["orb_diff_acr"].itrs[:, 0] ** 2 + dset["orb_diff_acr"].itrs[:, 1] ** 2))
    )
    dset.add_float("dradial", val=dset["orb_diff_acr"].itrs[:, 2])
    dset.add_float("clk_diff_dt_mean", val=dset["clk_diff"] - dset["clk_diff_with_dt_mean"])
    dset.write_as(tech=dset.vars["tech"], stage="where")

    return dset


def _generate_ure_control_tool_dset(inpath, rundate):
    """Generate URE control tool dataset

    Read URE control tool output file and generate afterwards Dataset.
    """
    log.info("Read URE control tool output file and generate afterwards Dataset.")
    p = parsers.parse_file(parser_name="ure_control_tool_csv", file_path=inpath)

    # Define dataset
    dset = data.Dataset(rundate=rundate, tech="sisre", stage="ure_tool", dataset_name="", dataset_id=0, empty=True)
    p.write_to_dataset(dset)

    # Write Dataset to file
    dset.write()

    return dset


############################################################################################
#
#
#                                MAIN PROGRAM
#
#
############################################################################################
if __name__ == "__main__":
    """
    """
    dset_diff_merged = None
    log.init(log_level="info")
    year = 2019
    month = 3
    day = 1
    rundate = date(year, month, day)
    nav_type = "FNAV"  # INAV, FNAV, LNAV
    system = "E"  # E, G
    frequencies = "E:E1_E5a"  # E:E1 (INAV), E:E1_E5b (INAV), E:E1_E5a (FNAV), G:L1_L2 (LNAV)
    if nav_type == "FNAV":
        inpath = "./data/G_GAL258_E1E5a_URE-AllPRN_{:%y%m%d}.csv".format(rundate)

    # Generate datasets
    dset_ure_tool = _generate_ure_control_tool_dset(inpath, rundate)
    dset_where = _generate_where_dset(rundate, nav_type, system, frequencies)

    # Generate difference between datasets
    _difference(dset_where, dset_ure_tool)
