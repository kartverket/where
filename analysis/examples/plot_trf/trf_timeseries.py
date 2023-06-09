#!/usr/bin/env python3
"""Plot station coordinates for one station in a set of reference frames

Usage:

    ./trf_timeseries.py tech station-id [year-from:year-to [trf_1 trf_2 ...]]

By default, the time series are plotted for the following reference frames:

    {trfs}

Examples:
---------

Plot default time period and list of trfs for VLBI-7345 (Tsukuba)

    ./trf_timeseries.py vlbi 7345

Plot VLBI-7345 from 2007 to 2014 (inclusive)

    ./trf_timeseries.py vlbi 7345 2007:2014

Compare ITRF2008 and ITRF2014 (SSC-files) for  VLBI-7345 in 2011

    ./trf_timeseries.py vlbi 7345 2011:2011 itrf:2008_ssc itrf:2014_ssc

Authors:
--------

* Geir Arne Hjelle <geir.arne.hjelle@kartverket.no>
"""

# Standard library imports
import sys

# External library imports
import numpy as np
import matplotlib.pyplot as plt

# Where imports
from where import apriori
from where.lib import config
from where.lib import log
from where.data.time import Time

# Reference frames that will be plotted
TRFS = ("itrf:2008_ssc", "itrf:2014_ssc", "itrf:2014_snx")


def read_trfs(station, time):
    """Read time series of positions for all TRFs for the given station"""
    pos = dict(time=time, trfs=list())

    # Read time series of positions from apriori TRFs
    for trf in TRFS:
        site = apriori.get("trf", time=time, reference_frames=trf)[station]
        trf_name = site.source
        pos["trfs"].append(trf_name)
        pos[trf_name] = site.pos.trs
        print("Reading {} from file".format(trf_name))

        # Remove missing positions (i.e. positions deep inside earth)
        missing_idx = np.linalg.norm(pos[trf_name], axis=1) < 5e6
        pos[trf_name][missing_idx] = np.nan

        # Normalize (all trfs are normalized by the same translation)
        if "norm" not in pos:
            pos["norm"] = np.nanmean(pos[trf_name], axis=0)
        pos[trf_name] -= pos["norm"]

    return pos


def plot_trfs(station, pos):
    """Plot X, Y, Z for all TRFs"""
    print("Plotting {} for {}".format(", ".join(pos["trfs"]), station.upper()))
    lines = dict()
    fig, axes = plt.subplots(3, sharex=True)
    for plot_num, xyz in enumerate("xyz"):
        axes[plot_num].get_yaxis().get_major_formatter().set_useOffset(False)
        axes[plot_num].get_yaxis().get_major_formatter().set_scientific(True)
        axes[plot_num].set_ylabel("{} [m]".format(xyz.upper()))
        for trf in pos["trfs"]:
            lines[trf] = axes[plot_num].plot(pos["time"].datetime, pos[trf][:, plot_num])[0]
    fig.legend([lines[t] for t in pos["trfs"]], pos["trfs"])
    fig.suptitle(station.upper())
    plt.show()


def _parse_time(time_str="1995:2018"):
    """Return a daily Time object for the given years (inclusive)"""
    pos_args, _ = _split_args()
    if len(pos_args) >= 3:
        time_str = pos_args[2]

    start_year, end_year = time_str.split(":")
    start = Time(float(start_year), fmt="decimalyear", scale="utc")
    end = Time(float(end_year) + 1, fmt="decimalyear", scale="utc")

    return Time(np.arange(start.mjd, end.mjd), fmt="mjd", scale="utc")


def _split_args(args=None):
    """Split command line arguments in positional and optional arguments"""
    if args is None:
        args = sys.argv[1:]

    pos_args, opt_args = list(), list()
    for arg in args:
        if arg.startswith("-"):
            opt_args.append(arg)
        else:
            pos_args.append(arg)

    return pos_args, opt_args


def main():
    log.init()
    pos_args, opt_args = _split_args()

    # TRFs given on the command line
    if len(pos_args) > 3:
        global TRFS
        TRFS = tuple(pos_args[3:])

    # Help
    if len(pos_args) < 2 or "-h" in opt_args or "--help" in opt_args:
        print(__doc__.format(trfs=", ".join(TRFS)))
        sys.exit()

    # Set technique (this should be improved ...)
    tech, station = pos_args[:2]
    config.set_analysis(rundate=None, tech=tech, session="")
    config.files.profiles = [tech]

    # Read and plot
    time = _parse_time()
    pos = read_trfs(station, time)
    plot_trfs(station, pos)


if __name__ == "__main__":
    main()
