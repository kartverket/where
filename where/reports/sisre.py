"""Write report about a SISRE analysis run (only for one session and not several sessions (stations))

Description:
------------

asdf




"""
# Standard library imports
from datetime import datetime
import re
import textwrap

# External library imports
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# WHERE imports
import where
from where import apriori
from where import cleaners
from where.lib import config
from where.lib import files
from where.lib import gnss
from where.lib import plugins
from where.reports import report
from where.writers import sisre_report


@plugins.register
def write_session_report(rundate, tech):

    dsets = dict()
    for _, dset, report_data in report.reports("remover_data"):
        remover_name = report_data["remover_name"]
        station = dset.dataset_name
        dset_station = dsets.setdefault("removers", dict()).setdefault(station, dset)
        dset_station.add_float(
            "keep_idx_{}".format(remover_name), val=report_data["keep_idx"]
        )  # TODO: What is keep_idx?

    for _, dset, report_data in report.reports("orbit_data"):
        dsets.setdefault("orbit", dict())[report_data["station"]] = dset

    session = list(dsets["orbit"].keys())[0]
    if not session:
        return 0

    dset = dsets["orbit"][session]
    sisre_report.sisre_report(dset)

    # with files.open("output_gnss_session_report", mode="wt") as fid:
    #     rejected_satellites_per_station(fid, dsets)
    #     rejected_satellite_observations(fid, dsets)
    #     rejected_station_observations(fid, dsets)


def rejected_satellites_per_station(fid, dsets):

    for dset in dsets.values():
        fid.write(
            "\n# Overview over rejected observations for all satellites for station {}\n\n"
            "".format(dset.dataset_name.upper())
        )
        fid.write("| Sat.  |   #obs |  #eobs | #eobs% |\n")
        fid.write("|-------|-------:|-------:|-------:|\n")
        for sat in dset.unique("satellite"):
            sat_idx = dset.filter(satellite=sat)
            num_obs = sum(sat_idx)
            num_edit_obs = sum(np.logical_not(dset.keep_idx[sat_idx]))
            fid.write(
                "| {sat:5s} | {obs:6d} | {eobs:6d} | {eobsp:6.0f} |\n"
                "".format(sat=sat, obs=num_obs, eobs=num_edit_obs, eobsp=100 * num_edit_obs / num_obs)
            )

        num_edit_obs = sum(np.logical_not(dset.keep_idx))
        fid.write(
            "| Sum   | {obs:6d} | {eobs:6d} | {eobsp:6.0f} |\n"
            "".format(obs=dset.num_obs, eobs=num_edit_obs, eobsp=100 * num_edit_obs / dset.num_obs)
        )


def rejected_satellite_observations(fid, dsets):

    sum_num_obs = sum_num_edit_obs = 0
    fid.write("\n# Overview over rejected observations for all satellites\n\n")
    fid.write("|   Sat.    |     #obs   |    #eobs   |   #eobs%   |\n")
    fid.write("|-----------|-----------:|-----------:|-----------:|\n")
    for dset in dsets.values():
        sum_num_obs += dset.num_obs
        num_edit_obs = sum(np.logical_not(dset.keep_idx))
        sum_num_edit_obs += num_edit_obs
        fid.write(
            "|   {sta:5s}   |   {obs:6d}   |   {eobs:6d}   |   {eobsp:6.0f}   |\n"
            "".format(
                sta=dset.dataset_name, obs=dset.num_obs, eobs=num_edit_obs, eobsp=100 * num_edit_obs / dset.num_obs
            )
        )
    fid.write(
        "| __Sum__   | __{obs:6d}__ | __{eobs:6d}__ | __{eobsp:6.0f}__ |\n"
        "".format(obs=sum_num_obs, eobs=sum_num_edit_obs, eobsp=100 * sum_num_edit_obs / sum_num_obs)
    )


def rejected_station_observations(fid, dsets):
    stations = sorted(dsets.keys())
    if not stations:
        return

    satellites = dsets[stations[0]].unique("satellite")  # TODO: Better?
    fid.write("\n# Overview over rejected observations for all stations\n\n")
    fid.write("| Sta.  |   #obs |  #eobs | #eobs% |\n")
    fid.write("|-------|-------:|-------:|-------:|\n")
    for sat in satellites:
        sum_num_obs = sum_num_edit_obs = 0
        for station in stations:
            dset = dsets[station]
            sat_idx = dset.filter(satellite=sat)
            num_obs = sum(sat_idx)
            sum_num_obs += num_obs
            num_edit_obs = sum(np.logical_not(dset.keep_idx[sat_idx]))
            sum_num_edit_obs += num_edit_obs
        fid.write(
            "| {sat:5s} | {obs:6d} | {eobs:6d} | {eobsp:6.0f} |\n"
            "".format(sat=sat, obs=sum_num_obs, eobs=sum_num_edit_obs, eobsp=100 * sum_num_edit_obs / sum_num_obs)
        )
