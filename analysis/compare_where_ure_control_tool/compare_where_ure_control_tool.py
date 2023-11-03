#!/usr/bin/env python3
"""Read URE control tool CSV output file and convert them to Where dataset

Description:
------------

"""
# Standard library imports
from datetime import date
from pathlib import Path

# External library imports
import numpy as np

# Midgard imports
from midgard import parsers
from midgard.data import dataset
from midgard.dev import log


def _difference(dset_where, dset_ure_tool):
    
    file_path = "./work/diff_where_ure_control_tool"

    if dset_where.num_obs == 0 or dset_ure_tool == 0:
        log.warn(
            f"Nothing to compare. Number of observations are zero at least for one dataset "
            f"(dset_where: {dset_where.num_obs}, dset_ure_tool: {dset_ure_tool.num_obs})."
        )
        return 1

    ddiff = dset_where.difference(dset_ure_tool, index_by="time, satellite")
    log.info(f"Write dataset with difference of Where and URE control tool fields: {file_path}")
    ddiff.write(file_path)


def _generate_where_dset(inpath, rundate):
    """Generate Where dataset

    """
    log.info(f"Read Where dataset file: {inpath}")
    
    if not inpath.exists():
        log.fatal(f"File path {inpath} does not exists.")
    
    # Define dataset
    dset = dataset.Dataset().read(inpath)
    dset.add_float(
        "sqrt_a2_c2", val=(np.sqrt(dset.orb_diff.acr.along ** 2 + dset.orb_diff.acr.cross ** 2))
    )
    dset.add_float("dradial", val=dset.orb_diff.acr.radial)
    dset.add_float("clk_diff_dt_mean", val=dset["clk_diff"] - dset["clk_diff_with_dt_mean"])
    dset.write("./result/gnss-{:%Y%m%d}-where-raw-dataset.hdf5".format(rundate))

    return dset


def _generate_ure_control_tool_dset(inpath, rundate):
    """Generate URE control tool dataset

    Read URE control tool output file and generate afterwards Dataset.
    """
    log.info(f"Read URE control tool output file: {inpath}")
    
    if not inpath.exists():
        log.fatal(f"File path {inpath} does not exists.")
        
    p = parsers.parse_file(parser_name="ure_control_tool_csv", file_path=inpath)
    dset = p.as_dataset()
    dset.write("./result/gnss-{:%Y%m%d}-uretool-raw-dataset.hdf5".format(rundate))

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
    year = 2023
    month = 5
    day = 1
    rundate = date(year, month, day)
    
    where_path = Path("./data/sisre-{:%Y%m%d}-calculate-raw-dataset.hdf5".format(rundate))
    #ure_path = Path("./data/G_GAL258_E1E5a_URE-AllPRN_{:%Y%m%d}.csv".format(rundate))
    ure_path = Path("./data/G_GPS_LNAVP1-P2_model_{:%Y%m%d}_120_None_None_OS.csv".format(rundate))

    # Generate datasets
    dset_ure_tool = _generate_ure_control_tool_dset(ure_path, rundate)
    dset_where = _generate_where_dset(where_path, rundate)

    # Generate difference between datasets
    _difference(dset_where, dset_ure_tool)
