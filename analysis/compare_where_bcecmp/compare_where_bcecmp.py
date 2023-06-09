#!/usr/bin/env python3
"""Compare Where and DLR BCEcmp SISRE analysis

Description:
------------
1. Read BCEcmp output file and generate afterwards Dataset.
2. Start SISRE analysis with Where and read Dataset.
3. Compare Where and BCEcmp Dataset

"""
# Standard library imports
from datetime import date
import subprocess

# External library imports
import numpy as np

# Where imports
from where import data
from midgard import parsers
from where.lib import log


def _difference_where_bcecmp(dwhr, dbce):

    log.info("Generate difference between SISRE Where and BCEcmp analysis.")
    if dwhr.num_obs == 0:
        log.fatal("Nothing to compare. Where dataset is empty.")
    if dbce.num_obs == 0:
        log.fatal("Nothing to compare. BCEcmp dataset is empty.")

    # Get common dataset
    hash_where = [(t.gps.isot + sn) for t, sn in dwhr.values("time", "satellite")]
    hash_bcecmp = [(t.gps.isot + sn) for t, sn in dbce.values("time", "satellite")]
    hash_both = sorted(set(hash_where) & set(hash_bcecmp))
    idx_w = [hash_where.index(h) for h in hash_both]
    idx_b = [hash_bcecmp.index(h) for h in hash_both]

    ddiff = data.Dataset(rundate=date, tech="sisre", stage="dwhr_bce", dataset_name="", dataset_id=0, empty=True)
    ddiff.num_obs = len(hash_both)

    common_fields = ["clk_diff_with_dt_mean", "orb_diff_3d", "sisre", "used_iode"]
    for field in common_fields:
        if field in dwhr.fields:
            if field in dbce.fields:
                diff = dwhr[field][idx_w] - dbce[field][idx_b]
                ddiff.add_float(field, val=diff)

    ddiff.add_float(
        "sqrt_a2_c2",
        val=(
            np.sqrt(dwhr["orb_diff_acr"].itrs[:, 0][idx_w] ** 2 + dwhr["orb_diff_acr"].itrs[:, 1][idx_w] ** 2)
            - np.sqrt(dbce["orb_diff_acr"].itrs[:, 0][idx_b] ** 2 + dbce["orb_diff_acr"].itrs[:, 1][idx_b] ** 2)
        ),
    )
    ddiff.add_time("time", val=dwhr.time.gps[idx_w], scale="gps")
    ddiff.add_position("orb_diff_acr", time="time", itrs=dwhr.orb_diff_acr.itrs[idx_w] - dbce.orb_diff_acr.itrs[idx_b])
    ddiff.add_text("satellite", val=dwhr.satellite[idx_w])
    ddiff.add_text("system", val=dwhr.satellite[idx_w].astype("U1"))

    ddiff.write()


def _generate_bcecmp_dset(inpath, date):
    """Generate BCEcmp dataset

    Read BCEcmp output file and generate afterwards Dataset.
    """
    log.info("Read BCEcmp output file and generate afterwards Dataset.")
    parser = parsers.parse_file(parser_name="bcecmp_sisre", file_path=inpath)
    bce_dict = parser.as_dict()

    # Define dataset
    dbce = data.Dataset(rundate=date, tech="sisre", stage="bcecmp", dataset_name="", dataset_id=0, empty=True)
    dbce.num_obs = len(bce_dict["time"])

    dbce.add_time("time", val=bce_dict["time"], scale="gps")
    dbce.add_text("satellite", val=bce_dict["satellite"])
    dbce.add_text("system", val=bce_dict["system"])

    for field in bce_dict.keys():
        if field not in ["time", "satellite", "system"]:
            dbce.add_float(field, val=bce_dict[field])

    dbce.add_position(
        "orb_diff_acr",
        time="time",
        itrs=np.vstack((dbce.dalong_track, dbce.dcross_track, dbce.dradial)).T,
        unit="meter",
    )
    dbce.add_float("orb_diff_3d", val=np.linalg.norm(dbce.orb_diff_acr.itrs, axis=1), unit="meter")

    # Write Dataset to file
    dbce.write()

    return dbce


def _generate_where_dset(date, nav_type, system, frequencies):
    """Generate SISRE dataset

    Start SISRE analysis with Where and read Dataset.
    """
    program = "where"
    tech = "sisre"
    session_name = "brdm"

    # Start processing with Where
    where_args = [
        str(date.year),
        str(date.month),
        str(date.day),
        "--" + tech,
        "--session_name=" + session_name,
        "--brdc_block_nearest_to=transmission_time:positive",
        "--systems=" + system,
        "--frequencies=" + frequencies,
        "--navigation_message_type=" + system + ":" + nav_type,
        "--clock_product=sp3",
        "--weight_factor_elev_mask=5",
        "--profile=bcecmp",
        # "-N",
        # "-D",
        "-T",
    ]
    log.info("Start {}: {} {}".format(tech.upper(), program, " ".join(where_args)))
    process = subprocess.run([program] + where_args)
    if process.returncode:
        log.error(f"Where {tech.upper()} failed with error code {process.returncode} ({' '.join(process.args)})")

    # Write SISRE Dataset
    dset_vars = dict(rundate=date, tech="sisre", stage="calculate", session="", dataset_name="", dataset_id="last")

    dset = data.Dataset(**dset_vars)
    dset.write_as(tech=tech, stage="where", dataset_name="")

    return dset


############################################################################################
#
#
#                                MAIN PROGRAM
#
#
############################################################################################
if __name__ == "__main__":

    log.init(log_level="info")
    # year = 2018
    year = 2019
    # month = 2
    # month = 3
    month = 8
    # day = 1
    # day = 31
    day = 28
    rundate = date(year, month, day)
    nav_type = "INAV"  # INAV, FNAV, LNAV
    system = "E"  # E, G
    frequencies = "E:E1"  # E:E1 (INAV), E:E1_E5b (INAV), E:E1_E5a (FNAV), G:L1_L2 (LNAV)
    if nav_type == "FNAV":
        inpath = "BCEcmp_GAL_{:s}_E1E5A_com_{:%Y}_{:%0j}.OUT".format(nav_type, rundate, rundate)
    elif nav_type == "INAV" and frequencies == "E:E1":
        # inpath = "BCEcmp_GAL_{:s}_E1_com_{:%Y}_{:%0j}.OUT".format(nav_type, rundate, rundate)
        inpath = "BCEcmp_GAL_{:s}_E1_COD0MGXFIN_{:%Y}_{:%0j}.OUT".format(nav_type, rundate, rundate)
    elif nav_type == "INAV" and frequencies == "E:E1_E5b":
        inpath = "BCEcmp_GAL_{:s}_E1E5B_com_{:%Y}_{:%0j}.OUT".format(nav_type, rundate, rundate)
    elif nav_type == "LNAV":
        inpath = "BCEcmp_GPS_{:s}_L1L2PY_com_{:%Y}_{:%0j}.OUT".format(nav_type, rundate, rundate)

    dbce = _generate_bcecmp_dset(inpath, rundate)
    dwhr = _generate_where_dset(rundate, nav_type, system, frequencies)
    _difference_where_bcecmp(dwhr, dbce)
