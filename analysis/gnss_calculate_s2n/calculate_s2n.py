#!/usr/bin/env python3
"""Calculate Signal-to-Noise based on a Rinex file

Usage:
------
    ./calculate_s2n.py rinex-files

Output is written to Where datasets (open with There), a CSV-file named s2n.csv and PNG-images in the current
directory.

Authors:
-------

* Geir Arne Hjelle <geir.arne.hjelle@kartverket.no>
"""

import itertools
import pathlib
import sys

import pandas as pd
import matplotlib.pyplot as plt

from where import data
from where.lib import log
from where import parsers


def main():
    log.init()

    # Read file paths from command line
    file_paths = [pathlib.Path(a) for a in sys.argv[1:] if not a.startswith("-")]
    if not file_paths:
        print(__doc__)
        sys.exit(1)

    # Collect data from rinex files and store in datasets
    s2n_data = list()
    for file_path in file_paths:
        s2n_data.extend(convert_rinex_to_dset(file_path))

    # Convert summary data to a dataframe and plot it
    s2n_df = pd.DataFrame(s2n_data, columns=s2n_data[0].keys())
    s2n_df.to_csv("s2n.csv")
    plot_df(s2n_df)


def convert_rinex_to_dset(file_path):
    log.info("Working with {}", file_path)

    # Read Rinex data and save them to the dataset
    dset = data.Dataset.anonymous()
    rnx_data = parsers.parse_file("rinex2_obs", file_path, sampling_rate=30)
    rnx_data.write_to_dataset(dset)

    # Add field for seconds
    dset.add_float("sec_of_day", val=dset.time.gps.sec_of_day, unit="seconds")

    # Find metadata for the dataset and write it to disk
    site = dset.station[0]
    date = dset.time.gps.datetime[0].date()
    tech = "gnss"
    stage = "s2n"
    dset.write_as(rundate=date, tech=tech, stage=stage, dataset_name=site)

    # Print and return some statistics about s2n
    return s2n_statistics(dset)


def s2n_statistics(dset):
    fields = ("S1", "S2")
    satellites = dset.unique("satellite", system="G")

    # Print statistics to screen
    print(f"\n{dset.dataset_name} {dset.rundate}\n" + "=" * 70)
    for satellite in satellites:
        print(f"{satellite:<3s}  " + "  ".join(f"{f} = {dset.mean(f, satellite=satellite):10.6f}" for f in fields))
    print("all  " + "  ".join(f"{f} = {dset.mean(f, system='G'):10.6f}" for f in fields))

    # Store statistics in a list of dicts that can be converted to Pandas Dataframes
    sat_data = [{s: dset.mean(f, satellite=s) for s in satellites} for f in fields]
    return [dict(station=dset.dataset_name, date=dset.rundate, obs=f, **s) for s, f in zip(sat_data, fields)]


def plot_df(df):
    for station, obs in itertools.product(df.station.unique(), df.obs.unique()):
        file_path = f"{station}_{obs}.png"
        log.info(f"Plotting {obs} for {station} to {file_path}")
        ax = df[(df.station == station) & (df.obs == obs)].plot(x="date", title=f"{station} - {obs}", legend=False)
        for label in ax.get_xticklabels():
            label.set_rotation(12)
        plt.savefig(file_path)


if __name__ == "__main__":
    main()
