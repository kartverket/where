"""Compare different SISRE Where datasets and write a text file in GRC CSV format

Description:
------------

A dictionary with datasets is used as input for this writer. The keys of the dictionary are station names. 

Example:
--------


"""
# Standard library imports
import csv
from typing import Any, Dict, List, Tuple

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.data import position
from midgard.dev import plugins

# Where imports
import where
from where.lib import config, log, util
from where.writers._gnss import get_grc_csv_header, get_grc_csv_row


@plugins.register
def sisre_comparison_grc_csv(dset: "Dataset") -> None:
    """Compare different SISRE Where datasets and write a text file in GRC CSV format

    Args:
        dset:  Dictionary with station name as keys and the belonging Dataset as value
    """
    output_defs = dict()
    dset_first = dset[list(dset.keys())[0]]
    _,nav_msg,mode = dset_first.vars["profile"].split("_")[0:3]

    # Get file path
    file_vars = {**dset_first.vars, **dset_first.analysis}
    file_vars["solution"] = config.tech.gnss_comparison_report.solution.str.lower()
    file_path = config.files.path("output_sisre_comparison_grc_csv", file_vars=file_vars)
    file_path.parent.mkdir(parents=True, exist_ok=True)

    # Get dataframes for writing
    df, df_month = _generate_dataframe(dset)

    # Write file
    # 
    # Example:
    #
    # Service Line,Service Category,Business Service,Batch,Satellite,PRN,Slot,GSS Site,Service,Type,Mode,Target,Unit,Month/Year,Result
    #   Open Service,Ranging Domain,SDD-OS-0012 SIS Ranging Accuracy over All Satellites over Month,,,,,,,OS,Single,E1,2,m,21-Jul,0.14685266875382683
    #   Open Service,Ranging Domain,SDD-OS-0012 SIS Ranging Accuracy over All Satellites over Month,,,,,,,OS,Dual,E1/E5b,2,m,21-Jul,0.12860206883303776
    #   Open Service,Ranging Domain,SDD-OS-0012 SIS Ranging Accuracy over All Satellites over Month,,,,,,,OS,Dual,E1/E5a,2,m,21-Jul,0.1318876874345714
    #   Open Service,Ranging Domain,SDD-OS-0013 SIS Ranging Accuracy per Satellite over Month,FOC,E01,GSAT210,A02,,,OS,Single,E1,7,m,21-Jul,0.22263083688024857
    #   Open Service,Ranging Domain,SDD-OS-0013 SIS Ranging Accuracy per Satellite over Month,FOC,E02,GSAT211,A06,,,OS,Single,E1,7,m,21-Jul,0.29602410930897965
    #
    with open(file_path, "w") as csvfile:
        writer = csv.writer(csvfile, delimiter=";")

        log.info(f"Write file {file_path}.")

        # Get signal combination modes
        signal_modes= list(df_month.columns)
        signal_modes.remove("date")

        # Write header 
        writer.writerow(get_grc_csv_header())

        # Write constellation results to GRC CSV file
        for mode in signal_modes:

            # Get constellation name
            if mode.startswith("E"):
                constellation = "Galileo"
            elif mode.startswith("L"):
                constellation = "GPS"

            for _, row in df_month.iterrows():

                writer.writerow(
                    get_grc_csv_row(
                        constellation=constellation,
                        kpi="sisre", 
                        mode=mode, 
                        date=row.date, 
                        result=row[mode],
                    )
                )

        # Write satellite-wise results to GRC CSV file
        for mode in signal_modes:

            # Filtering necessary to get only satellites related to one GNSS
            if mode.startswith("E"):
                system = "E"
                constellation = "Galileo"
            elif mode.startswith("L"):
                system = "G"
                constellation = "GPS"

            idx = df.system == system

            # Get satellites
            satellites = set(df[idx].satellite)

            df_mode = df[idx].pivot(index="time_gps", columns="satellite", values=mode)
            df_mode = df_mode.resample("M").apply(lambda x: np.nanpercentile(x, q=95))
            df_mode.index = df_mode.index.strftime("%Y-%b")
            df_mode.index.name = "date"
            df_mode = df_mode.reset_index()

            for index, row in df_mode.iterrows():
                for sat in sorted(satellites):
                    writer.writerow(
                        get_grc_csv_row(
                            constellation=constellation,
                            kpi="sisre_sat", 
                            mode=mode, 
                            date=row.date, 
                            result=row[sat],
                            satellite=sat,
                        )
                    )
                


def _generate_dataframe(dsets: Dict[str, "Dataset"]) -> Tuple[pd.core.frame.DataFrame]:
    """Generate dataframes based on SISRE datasets

    The dataframe "df" has following columns:

        time_gps:       Time in GPS time scale given as datetime objects
        satellite:      Satellite identifiers
        system:         GNSS identifier
        <solution_1>:   First SISRE solution (e.g. E1)
        <solution_2>:   Second SISRE solution (e.g. E1/E5b)
        <solution_3>:   Second SISRE solution (e.g. E1/E5a)

    Example for "df_month_perc_rms" dictionary:

                        E1    E1/E5b    E1/E5a
        2021-Jul  0.335688  0.297593  0.326859
        2021-Jul  0.380575  0.330701  0.352535
        2021-Jul  0.353586  0.314817  0.344597


    Args:
        dsets: Dictionary with SISRE solution name as keys (e.g. cnes_inav_e1, cnes_inav_e1e5b, cnes_fnav_e1e5a) and
               the belonging Dataset as value

    Returns:

        Tuple with following entries:

        | Element              | Description                                                                          |
        |----------------------|--------------------------------------------------------------------------------------|
        | df                   | Given DAILY SISRE solutions are merged into one dataframe                            |
        | df_month_perc_rms    | Dataframe with MONTHLY samples of 95th percentile SISRE, which are based on epochwise|
        |                      | RMS SISRE solutions (based on Galileo SDD v1.1 version)                              |

    """
    df = pd.DataFrame()
    df_month_perc_rms = pd.DataFrame()
    signal_types = list()

    for name, dset in dsets.items():

        if dset.num_obs == 0:
            log.warn(f"Dataset '{name}' is empty.")
            continue

        signal_type = _get_signal_type(dset.meta)
        signal_types.append(signal_type)
        df_tmp = dset.as_dataframe(fields=["satellite", "system", "sisre", "time.gps"])  # , index="time.gps")
        df_tmp = df_tmp.rename(columns={"sisre": signal_type})

        if df.empty:
            df = df_tmp
            continue
        df = df.merge(df_tmp, on=["satellite", "system", "time_gps"], how="outer")

    if df.empty:
        log.fatal(f"All given datasets are empty [{', '.join(dsets.keys())}].")

    # Generate monthly samples of 95th percentile SISRE based on epochwise SISRE RMS solutions(after SDD v1.1 version)
    epochs = sorted(set(df["time_gps"]))
    df_tmp = pd.DataFrame(index=epochs, columns=signal_types)

    # Loop over observation epochs
    for epoch in epochs:
        idx = df["time_gps"] == epoch
        row = dict()

        # Determine RMS for each signal type over all given SISRE satellite solutions in each epoch
        for signal_type in signal_types:
            row[signal_type] = np.sqrt(np.nanmean(np.square(df[signal_type][idx])))
        df_tmp.loc[epoch] = pd.Series(row)

    df_month_perc_rms = df_tmp.resample("M").apply(lambda x: np.nanpercentile(list(x), q=95))
    df_month_perc_rms.index = df_month_perc_rms.index.strftime("%Y-%b")

    # Prepare dataframes for writing
    df_month_perc_rms.index.name = "date"
    df_month_perc_rms = df_month_perc_rms.reset_index()

    return df, df_month_perc_rms


def _get_signal_type(meta: Dict[str, Any]) -> str:
    """Get signal type used for SISRE dataset

    Args:
        meta:   Dataset meta dictionary

    Returns:
        Signal type (e.g. E1, E1/E5a, ...)
    """
    if len(meta["systems"]) > 1:
        log.fatal(
            f"The writer '{FILE_NAME}' can only be used, if the dataset is based on one GNSS "
            f"(not '{', '.join(meta['systems'])})'."
        )

    system = meta["systems"][0]
    try:
        signal_type = meta["frequencies"][system].replace("_", "/")

    except KeyError:
        log.fatal(f"No frequencies are defined for GNSS {system!r} for option 'frequencies'.")

    return signal_type

