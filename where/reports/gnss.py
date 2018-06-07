"""Write report about a GNSS model run

Description:
------------

asdf




"""
# Standard library imports
from datetime import datetime

# External library imports
import numpy as np

# WHERE imports
import where
from where.lib import config
from where.lib import files
from where.lib import plugins
from where.reports import report


@plugins.register
def write_session_report(rundate, tech):
    dsets = dict()
    for _, dset, report_data in report.reports("remover_data"):
        dset.add_float("keep_idx", val=report_data["keep_idx"])
        dsets[dset.dataset_name] = dset

    with files.open("output_gnss_session_report", mode="wt") as fid:
        header(fid)
        rejected_satellites_per_station(fid, dsets)
        rejected_satellites(fid, dsets)
        rejected_stations(fid, dsets)


def header(fid):
    _, _, title = next(report.reports("title"))

    fid.write(
        """---
title: {title}
author: Where v{version} [{user}]
date: {nowdate:%Y-%m-%d}
---
""".format(
            title=title["text"], version=where.__version__, user=config.analysis.user.str, nowdate=datetime.now()
        )
    )

    fid.write("# bla bla bla\n")


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


def rejected_stations(fid, dsets):

    sum_num_obs = sum_num_edit_obs = 0
    fid.write("\n# Overview over rejected observations for all satellites\n\n")
    fid.write("|   Sta.    |     #obs   |    #eobs   |   #eobs%   |\n")
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
        "| __Sum__   | {obs:>10} | {eobs:>10} | {eobsp:>10} |\n"
        "".format(
            obs="__{}__".format(sum_num_obs),
            eobs="__{}__".format(sum_num_edit_obs),
            eobsp="__{:.0f}__".format(100 * sum_num_edit_obs / sum_num_obs),
        )
    )


def rejected_satellites(fid, dsets):
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
