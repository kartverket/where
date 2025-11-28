"""Compare different GNSS Where datasets and write a text file in GRC CSV format

Description:
------------

A dictionary with datasets is used as input for this writer. The keys of the dictionary are station names. 

Example:
--------


"""
# Standard library imports
import csv
from typing import Dict, List, Tuple

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.collections import enums
from midgard.data import position
from midgard.dev import plugins

# Where imports
import where
from where.lib import config, log
from where.writers._gnss import get_grc_csv_header, get_grc_csv_row


@plugins.register
def gnss_comparison_grc_csv(dset: "Dataset") -> None:
    """Compare different GNSS Where datasets and write a text file in GRC CSV format

    Args:
        dset:  Dictionary with station name as keys and the belonging Dataset as value
    """
    output_defs = dict()
    dset_first = dset[list(dset.keys())[0]]
    _,nav_msg,mode = dset_first.vars["profile"].split("_")[0:3]

    # Get file path
    file_vars = {**dset_first.vars, **dset_first.analysis}
    file_vars["solution"] = config.tech.gnss_comparison_report.solution.str.lower()
    file_path = config.files.path("output_gnss_comparison_grc_csv", file_vars=file_vars)
    file_path.parent.mkdir(parents=True, exist_ok=True)

    # Get dataframes for writing
    constellation, df_month = _generate_dataframe(dset)

    # Write file
    # 
    # Example:
    #
    #   Service Line,Service Category,Business Service,Batch,Satellite,PRN,Slot,GSS Site,Service,Type,Mode,Target,Unit,Month/Year,Result
    #   Open Service,Position Domain,FOM-OS-0015 Horizontal Positioning Service Accuracy per Station over Month,,,,,Ny Alesund (Norway),nabd,OS,Single,E1,7.5,m,July-2021,0.92873542411617
    #   Open Service,Position Domain,FOM-OS-0016 Vertical Positioning Service Accuracy per Station over Month,,,,,Ny Alesund (Norway),nabd,OS,Single,E1,15,m,July-2021,2.8155215806884515
    #
    with open(file_path, "w") as csvfile:
        writer = csv.writer(csvfile, delimiter=";")

        log.info(f"Write file {file_path}.")

        # Write header 
        writer.writerow(get_grc_csv_header())

        # Write results to GRC CSV file
        for kpi in ["hpe", "vpe"]:

            for index, row in df_month.iterrows():
                writer.writerow(
                    get_grc_csv_row(
                        constellation=constellation,
                        kpi=kpi, 
                        mode=mode, 
                        date=row.date, 
                        result=row[kpi],
                        station=row.station, 
                    )
                )


def _generate_dataframe(dset: Dict[str, "Dataset"]) -> Tuple[str, Dict[str, pd.core.frame.DataFrame]]:
    """Generate dataframe based on station datasets

    Example for "df_month" dictionary:

               date           hpe           vpe station
        0  2021-Jul  1.936327e-09  5.567021e-09    nabd
        1  2021-Jul  3.034691e-09  5.419228e-09    hons
        2  2021-Jul  3.865243e-09  5.229739e-09    vegs

    Args:
        dset: Dictionary with station name as keys and the belonging Dataset as value

    Returns:
        Tuple with GNSS constellation information and dataframe with monthly entries with columns like date, station,
        hpe, vpe, ...

    """
    dsets = dset
    df_month = pd.DataFrame()
    dfs_month_field = {
        "hpe": pd.DataFrame(), 
        "vpe": pd.DataFrame(), 
    }

    for station, dset in dsets.items():

        if dset.num_obs == 0:
            log.warn(f"Dataset '{station}' is empty.")
            continue

        if dset.unique("system").size > 1:
            log.fatal(f"Data from more than one GNSS is given in dataset. The writer '{__name__}' can only "
                      f"handle data from one GNSS.")
        else:
            constellation = enums.gnss_id_to_name[dset.unique("system")[0]]
            
        # Determine topocentric coordinates (east, north, up)
        ref_pos = position.Position(
                            np.repeat(
                                  np.array([dset.meta["pos_x"], dset.meta["pos_y"], dset.meta["pos_z"]])[None, :], 
                                  dset.num_obs, 
                                  axis=0,
                            ), 
                            system="trs",
        )

        if not "enu" in dset.fields:
            dset.add_position_delta(
                name="enu",
                val=(dset.site_pos.trs - ref_pos).val,
                system="trs",
                ref_pos=ref_pos,
                write_level="operational",
            )

        # TODO: Maybe it is not necessary to introduce enu, hpe and vpe to dataset
        #      Maybe better to introduce fields in estimate stage already.
        if not "hpe" in dset.fields:
            hpe = np.sqrt(dset.enu.enu.east ** 2 + dset.enu.enu.north ** 2)
            dset.add_float("hpe", val=hpe, write_level="operational")

        if not "vpe" in dset.fields:
            vpe = np.absolute(dset.enu.enu.up)
            dset.add_float("vpe", val=vpe, write_level="operational")

        # Determine dataframe 
        df = dset.as_dataframe(fields=["time.gps", "hpe", "vpe"])
        df = df.rename(columns={"time_gps": "date"})

        if df.empty:
            continue
        else:
            df_month_tmp = df.set_index("date").resample("M").apply(lambda x: np.nanpercentile(x, q=95))
            df_month_tmp.index = df_month_tmp.index.strftime("%Y-%b")
            df_month_tmp["station"] = np.repeat(station, df_month_tmp.shape[0])
            df_month = pd.concat([df_month, df_month_tmp], axis=0)


    # Prepare dataframes for writing
    df_month.index.name = "date"
    df_month = df_month.reset_index()

    return constellation, df_month



