import argparse
from datetime import date, datetime, timedelta
import sys
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import matplotlib.ticker as mt
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
from scipy.spatial import Delaunay

from midgard.math.unit import Unit
from midgard.math.constant import constant

from where.data import dataset3 as dataset
from where.data.position import PositionDelta
from where.data.time import Time
from where.lib import config
from where.lib import rotation
from where import apriori


def plot(x, ys, errors, colors, labels, name, station):
    num_plots = len(ys)
    fig, axs = plt.subplots(num_plots, figsize=(12, 6), dpi=150, sharex=True)
    for i, (y, e, l) in enumerate(zip(ys, errors, labels)):
        axs[i].errorbar(x, y, yerr=e, fmt="o", marker=None, zorder=0, mew=0, ecolor="tab:gray")
        im = axs[i].scatter(x, y, c=colors, zorder=100)
        axs[i].set_ylim((-0.1, 0.1))
        axs[i].set(ylabel=l)
        axs[i].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        axs[i].xaxis.set_major_locator(mt.LinearLocator(numticks=7))
        axs[i].grid(axis="y", linestyle="--")
    cbar = fig.colorbar(im, ax=axs.ravel().tolist(), use_gridspec=True)
    cbar.set_label("Number of observations used in solution")
    fig.autofmt_xdate()
    #fig.suptitle(station)
    #plt.tight_layout()
    plt.savefig(f"img/{dset_id}/Position_{station}_{name}_{dset_id}.png", bbox_inches='tight')
    plt.close()


def plot_sta_param(x, ys, num, labels, title, ylabel):
    num_plots = len(ys)
    fig = plt.figure(figsize=(12, 6), dpi=150)
    xlim = (x.min(), x.max())
    for i, (y, n, l) in enumerate(zip(ys, num, labels)):
        legend_text = f"{l} (num_obs: {n})"
        plt.plot(x, y, marker=".", linestyle=None,label=legend_text)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    plt.gca().xaxis.set_major_locator(mt.LinearLocator(numticks=7))

    #plt.xlim(xlim)
    plt.legend()
    plt.ylabel(ylabel)
    plt.grid(axis="y", linestyle="--")
    #fig.suptitle(title)
    fig.autofmt_xdate()
    title = title.replace(" ", "_")
    sub_dir = "_".join(title.split("_")[0:-2])
    os.makedirs(f"img/{dset_id}/{sub_dir}", exist_ok=True)
    plt.savefig(f"img/{dset_id}/{sub_dir}/{title}_{dset_id}.png", bbox_inches='tight')
    plt.close()


def plot_skycoverage(dset_session, station):
    idx0 = dset_session.filter(station=station)
    idx1 = dset_session.filter(station_1=station)
    idx2 = dset_session.filter(station_2=station)
    elevation1 = dset_session.site_pos_1.elevation[idx1]
    elevation2 = dset_session.site_pos_2.elevation[idx2]
    azimuth1 = dset_session.site_pos_1.azimuth[idx1]
    azimuth2 = dset_session.site_pos_2.azimuth[idx2]
    residuals = dset_session.residual[idx0]
    elevation = np.concatenate((elevation1, elevation2))
    azimuth = np.concatenate((azimuth1, azimuth2))

    title = f"{station} {dset_session.vars['rundate']} {dset_session.vars['session_code']}"

    fig = plt.figure(figsize=(8, 8), dpi=150)
    ax = fig.add_subplot(111, projection="polar")
    im = ax.scatter(azimuth, elevation, c=residuals)
    ax.set_title(title)
    cbar = plt.colorbar(im, use_gridspec=True)
    cbar.set_label("Residual after estimation [m]")
    plt.tight_layout()
    title = title.replace(" ", "_")
    sub_dir = "Skyplot"
    os.makedirs(f"img/{dset_id}/{sub_dir}", exist_ok=True)
    fig.savefig(f"img/{dset_id}/{sub_dir}/{sub_dir}_{title}_{dset_id}.png")
    plt.close()


def plot_station_pos(dset_ts, station):

    # Select data from dataset
    idx = dset_ts.filter(station=station)
    dates = np.array([datetime.strptime(dt, "%Y-%m-%d") for dt in dset_ts.rundate[idx]])
    num_obs = np.zeros(np.sum(dset_ts.filter(station="all")))
    colors = dset_ts.num_obs_estimate[idx]
    all_idx = dset_ts.filter(station="all")
    all_dates = np.array([datetime.strptime(dt, "%Y-%m-%d") for dt in dset_ts.rundate[all_idx]])
    _, idx2, _ = np.intersect1d(all_dates, dates, return_indices=True)
    num_obs[idx2] = colors
    session_code = dset_ts.session_code[idx]

    t = Time(val=dates, scale="utc", fmt="datetime")
    trf = apriori.get("trf", time=t, reference_frames="itrf:2014, vtrf, custom")
    names = apriori.get("vlbi_station_codes")

    pos = trf[names[station]["cdp"]].pos
    lat, lon, _ = pos.llh.T
    trs2enu = rotation.enu2trs(lat, lon)
    enu2trs = rotation.trs2enu(lat, lon)

    dpos = PositionDelta(val=np.squeeze(dset_ts.neq_vlbi_site_pos[idx]), system="trs", ref_pos=pos)
    dpos_cov_xyz = dset_ts.neq_vlbi_site_pos_cov_[idx]
    dpos_ferr_xyz = np.sqrt(dpos_cov_xyz.diagonal(axis1=1, axis2=2))
    dpos_cov_enu = trs2enu @ dpos_cov_xyz @ enu2trs
    dpos_ferr_enu = np.sqrt(dpos_cov_enu.diagonal(axis1=1, axis2=2))

    # Plot timeseries of estimated station coordinate corrections
    plot(dates, dpos.val.T, dpos_ferr_xyz.T, colors, ["X [m]", " Y [m]", "Z [m]"], "xyz", station)
    plot(dates, dpos.enu.val.T, dpos_ferr_enu.T, colors, ["E [m]", " N [m]", "U [m]"], "enu", station)
    return num_obs, session_code


def plot_residuals(dset, sta_1, sta_2):
    plt.figure(figsize=(12, 8), dpi=150)
    plt.scatter(dset.time.utc.datetime, dset.residual, s=10)
    xlim = (min(dset.time.utc.datetime), max(dset.time.utc.datetime))
    plt.xlim(xlim)
    plt.ylabel("Postfit residuals [m]")
    plt.title(f"{dset.vars['rundate']} {dset.vars['session_code']}, Root Mean Square: {rms:6.4f} [m]")
    plt.gca().get_xaxis().set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    plt.gca().xaxis.set_major_locator(mt.LinearLocator(numticks=7))
    plt.gca().yaxis.get_major_formatter().set_useOffset(False)
    plt.grid(axis="y", linestyle="--")
    plt.tight_layout()
    sub_dir = "Postfit_Residuals"
    os.makedirs(f"img/{dset_id}/{sub_dir}", exist_ok=True)
    plt.savefig(f"img/{dset_id}/{sub_dir}/{sub_dir}_{dset.vars['rundate']}_{dset.vars['session_code']}_{dset_id}.png")
    plt.close()

    idx1 = dset.filter(station=sta_1)
    if any(idx1):
        plt.figure(figsize=(12, 8), dpi=150)
        plt.xlim(xlim)
        plt.ylim((-0.06, 0.06))
        plt.scatter(dset.time.utc.datetime[idx1], dset.residual[idx1])
        plt.title(f"{sta_1} {dset.vars['rundate']} {dset.vars['session_code']}, Root Mean Square: {rms1:6.4f} [m]")
        plt.ylabel("Postfit residuals [m]")
        plt.gca().get_xaxis().set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        plt.gca().xaxis.set_major_locator(mt.LinearLocator(numticks=7))
        plt.gca().yaxis.get_major_formatter().set_useOffset(False)
        plt.grid(axis="y", linestyle="--")
        plt.savefig(
            f"img/{dset_id}/{sub_dir}/{sub_dir}_{sta_1}_{dset.vars['rundate']}_{dset.vars['session_code']}_{dset_id}.png"
        )
        plt.close()

    idx2 = dset.filter(station=sta_2)
    if any(idx2):
        rms2 = dset.rms("residual", station=sta_2)
        plt.figure(figsize=(12, 8), dpi=150)
        plt.xlim(xlim)
        plt.ylim((-0.06, 0.06))
        plt.scatter(dset.time.utc.datetime[idx2], dset.residual[idx2])
        plt.title(f"{sta_2} {dset.vars['rundate']} {dset.vars['session_code']}, Root Mean Square: {rms2:6.4f} [m]")
        plt.ylabel("Postfit residuals [m]")
        plt.gca().get_xaxis().set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        plt.gca().xaxis.set_major_locator(mt.LinearLocator(numticks=7))
        plt.gca().yaxis.get_major_formatter().set_useOffset(False)
        plt.grid(axis="y", linestyle="--")
        plt.savefig(
            f"img/{dset_id}/{sub_dir}/Postfit_Residuals_{sta_2}_{dset.vars['rundate']}_{dset.vars['session_code']}_{dset_id}.png"
        )
        plt.close()


def plot_residual_rms(dates, rms, rms1, rms2, sta_1, sta_2, colors):
    num_plots = 3
    ys = [rms, rms1, rms2]
    labels = ["All stations", sta_1, sta_2]
    norm = mcolors.Normalize(np.min(np.concatenate(colors)), np.max(np.concatenate(colors)))
    im = cm.ScalarMappable(norm=norm)
    fig, axs = plt.subplots(num_plots, figsize=(12, 6), dpi=150, sharex=True)
    for i, (y, l) in enumerate(zip(ys, labels)):
        color = np.full(len(y), fill_value=np.nan)
        color[~np.isnan(y)] = colors[i][colors[i]!=0]
        axs[i].scatter(dates, y, c=color, norm=norm, zorder=100)
        axs[i].set(ylabel=l)
        axs[i].set_ylim((0.0, 0.03))
        #axs[i].set_xlim((min(dates) - timedelta(days=1), max(dates) + timedelta(days=1)))
        axs[i].grid(axis="y", linestyle="--")
        axs[i].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        axs[i].xaxis.set_major_locator(mt.LinearLocator(numticks=7))
        axs[i].grid(axis="y", linestyle="--")

    cbar = fig.colorbar(im, ax=axs.ravel().tolist(), use_gridspec=True)
    cbar.set_label("Number of observations used in solution")
    plt.suptitle("Root Mean Square of Postfit Residuals  [m]")
    #plt.tight_layout()
    fig.autofmt_xdate()
    plt.savefig(f"img/{dset_id}/RMS_Postfit_Residuals_{sta_1}_{sta_2}_{dset_id}.png", bbox_inches='tight')
    plt.close()


def plot_statistics(dates, dof, variance_factor, colors):
    fig, axs = plt.subplots(2, figsize=(12, 6), dpi=150, sharex=True)

    axs[0].set_ylabel('Degrees of freedom')
    #axs[0].set_ylim((1500, 14500))
    #axs[0].set_xlim((min(dates) - timedelta(days=1), max(dates) + timedelta(days=1)))
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    axs[0].xaxis.set_major_locator(mt.LinearLocator(numticks=7))
    im = axs[0].scatter(dates, dof, c=colors)
    axs[0].grid(axis="y", linestyle="--")
    
    axs[1].set_ylabel('Variance factor')
    #axs[1].set_ylim((0.5, 2.0))
    #axs[1].set_xlim((min(dates) - timedelta(days=1), max(dates) + timedelta(days=1)))
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    axs[1].xaxis.set_major_locator(mt.LinearLocator(numticks=7))
    axs[1].scatter(dates, variance_factor, c=colors)
    axs[1].grid(axis="y", linestyle="--")
    
    cbar = fig.colorbar(im, ax=axs.ravel().tolist(), use_gridspec=True)
    cbar.set_label("Number of observations used in solution")
    fig.autofmt_xdate()
    #fig.tight_layout()
    plt.savefig(f"img/{dset_id}/Statistics_{dset_id}.png", bbox_inches='tight')

def plot_baseline(dates, baseline_length, baseline_length_ferr, num_obs_bs, sta_1, sta_2, local_tie):
    t = Time(dates, fmt="datetime", scale="utc")
    _wblr, wmean = wblr(baseline_length,t.mjd, baseline_length_ferr)
    if local_tie:
        print(f"{dset_id}: Weighted mean: {wmean: 6.4f} [m], Local tie offset: {(wmean - local_tie)*Unit.m2mm: 6.2f} [mm]")
    else:
        print(f"{dset_id}: Weighted mean: {wmean: 16.4f} [m]")
    print(f"{dset_id}: Weighted blr (using weighted mean): {_wblr*Unit.m2mm: 6.2f} [mm]")

    x_padding = 1.5 # days
    y_padding = 0.01 # meter

    fig = plt.figure(figsize=(12, 8), dpi=150)
    xlim = (min(dates) - timedelta(days=x_padding), max(dates) + timedelta(days=x_padding))
    ylim = (wmean - y_padding, wmean + y_padding)
    plt.errorbar(dates, baseline_length, yerr=baseline_length_ferr, fmt="o", marker=None, zorder=0, mew=0, ecolor="tab:gray")
    im = plt.scatter(dates, baseline_length, c=num_obs_bs, zorder=100)
    plt.axhline(wmean, c="red",label="Weighted mean")
    if local_tie:
       plt.axhline(local_tie, c="blue", label="Local tie vector")
    #plt.xlim(xlim)
    plt.ylim(ylim)
    plt.title(f"Weighted Baseline Length Repeatability: {_wblr*Unit.meter2millimeter:4.2f} [mm]")
    plt.ylabel(f"Baseline length [m]")
    plt.legend()
    plt.gca().get_xaxis().set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    plt.gca().xaxis.set_major_locator(mt.LinearLocator(numticks=7))
    plt.gca().yaxis.get_major_formatter().set_useOffset(False)
    plt.grid(axis="y", linestyle="--")
    cbar = fig.colorbar(im, use_gridspec=True)
    cbar.set_label("Number of observations on baseline used in solution")
    #plt.tight_layout()
    fig.autofmt_xdate()
    plt.savefig(f"img/{dset_id}/Baseline_{sta_2}_{sta_1}_{dset_id}.png", bbox_inches='tight')
    plt.close()


def wblr(bl, mjd, ferr):
    """Computes weighted baseline length repeatability
    bl: estimated baseline lengths for one baseline [m]
    mjd: modified julian date for estimated baseline length
    w: formal error of estimated baseline lengths [m]

    Returns: weighted baseline length repeatability
    """
    bl = np.array(bl)
    mjd = np.array(mjd)
    ferr = np.array(ferr)

    b = bl
    w = 1 / ferr ** 2
    w_sum = np.sum(w)

    b_wmean = np.sum(w * b) / w_sum
    wblr1 = np.sqrt(np.sum((b - b_wmean) ** 2 * w) / (w_sum - np.sum(w ** 2) / w_sum))

    return wblr1, b_wmean


def plot_session_hist(num_obs_sta1, num_obs_sta2):
    
    def make_hist(data, station):
        fig = plt.figure(figsize=(12, 8), dpi=150)
        plt.bar(data["date"], data["schedule"], label="Scheduled", width=3)
        plt.bar(data["date"], data["read"], label="Recovered", width=3)
        plt.bar(data["date"], data["estimate"], label="Used", width=3)
        plt.legend()
        xlim = (data['date'].min()-timedelta(days=1), data['date'].max()+timedelta(days=1))
        plt.xlim(xlim)
        plt.gca().get_xaxis().set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        plt.gca().xaxis.set_major_locator(mt.LinearLocator(numticks=7))
        plt.gca().yaxis.get_major_formatter().set_useOffset(False)
        plt.ylabel("Number of observations")
        plt.xlim((data['date'].min(),data['date'].max()))
        fig.autofmt_xdate()
        plt.savefig(f"img/{dset_id}/Num_obs_{station}_{dset_id}.png", bbox_inches='tight')
        plt.close()

        with open(f"session_stats_{station}.txt", "w") as fid:
            fid.write(f"session code, date, scheduled, recovered, used\n")
            schedule = np.nan_to_num(data["schedule"])
            read = np.nan_to_num(data["read"])
            estimate = np.nan_to_num(data["estimate"])
            for sc, d, s, r, u in zip(data["session_code"], data["date"], schedule, read, estimate): 
                fid.write(f"{sc:8}, {d:%Y-%m-%d}, {int(s)}, {int(r)}, {int(u)} \n")

    make_hist(num_obs_sta1, station1)
    make_hist(num_obs_sta2, station2)

def get_state_values(dset, fieldname, fill_value=np.nan):
    try:
        values = dset.state[fieldname]
    except AttributeError:
        values = np.full(shape=(dset.num_obs), fill_value=fill_value)
    return values

# Program starts execution here
    

parser = argparse.ArgumentParser(epilog="Example: python plot_nyale13s.py --id nyale13s0 --stations NYALE13S NYALES20 --plot_skycoverage --plot_num_obs --plot_trop --plot_residuals")
parser.add_argument("--id", help="Dataset id of result files.", type=str, default="nyale13s0")
parser.add_argument("--stations", help="Name of the two stations in the baseline", nargs=2, type=str, default=["NYALE13S", "NYALES20"])
parser.add_argument("--plot_skycoverage", help="Enable this flag to plot sky coverage.", action="store_true")
parser.add_argument("--plot_num_obs", help="Enable this flag to plot num obs histogram.", action="store_true")
parser.add_argument("--plot_trop", help="Enable this flag to plot troposphere parameters for each session for specified stations", action="store_true")
parser.add_argument("--plot_residuals", help="Enable this flag to plot residuals for each session for specified stations", action="store_true")
args = parser.parse_args()

dset_id = args.id
station1, station2 = sorted(args.stations)

baseline1 = f"{station1}/{station2}"
baseline2 = f"{station2}/{station1}"
pipeline = "vlbi"

os.makedirs(f"img/{dset_id}", exist_ok=True)

# Setup config to make apriori modules work
config.set_analysis(rundate=None, pipeline=pipeline)
config.files.profiles = [pipeline]
config.read_pipeline(pipeline)

# Get data from timeseries dataset
dset_ts = dataset.Dataset.read(
    rundate=date(1970, 1, 1), pipeline=pipeline, stage="timeseries", label="0", session_code="", id=dset_id
)

idx = dset_ts.filter(station="all")
dates = np.array([datetime.strptime(dt, "%Y-%m-%d") for dt in dset_ts.rundate[idx]])
session_codes = dset_ts.session_code[idx]

colors = [dset_ts.num_obs_estimate[idx], 
          dset_ts.num_obs_estimate[dset_ts.filter(station=station1)],
          dset_ts.num_obs_estimate[dset_ts.filter(station=station2)]]



num_obs_bs = []
num_obs_sta1 = {}
num_obs_sta2 = {}
rms = []
rms_sta_1 = []
rms_sta_2 = []
baseline_length = []
baseline_length_ferr = []
baseline_dates = []

dof = []
variance_factor = []

idx_sta1 = dset_ts.filter(station=station1)
idx_sta2 = dset_ts.filter(station=station2)

num_obs_sta1["date"] = np.array([datetime.strptime(dt, "%Y-%m-%d") for dt in dset_ts.rundate[idx_sta1]])
num_obs_sta2["date"] = np.array([datetime.strptime(dt, "%Y-%m-%d") for dt in dset_ts.rundate[idx_sta2]])
num_obs_sta1["session_code"] = dset_ts.session_code[idx_sta1]
num_obs_sta2["session_code"] = dset_ts.session_code[idx_sta2]
num_obs_sta1["schedule"] = dset_ts.num_obs_schedule[idx_sta1]
num_obs_sta2["schedule"] = dset_ts.num_obs_schedule[idx_sta2]
num_obs_sta1["read"] = dset_ts.num_obs_read[idx_sta1]
num_obs_sta2["read"] = dset_ts.num_obs_read[idx_sta2]
num_obs_sta1["estimate"] = dset_ts.num_obs_estimate[idx_sta1]
num_obs_sta2["estimate"] = dset_ts.num_obs_estimate[idx_sta2]


# Plot per session
for rundate, session_code in zip(dates, session_codes):
    print(f"{rundate:%Y-%m-%d} {session_code}")

    dset_session = dataset.Dataset.read(
        rundate=rundate, pipeline=pipeline, stage="postprocess", label="last", session_code=session_code, id=dset_id
    )
    n1 = dset_session.num(station=station1)
    n2 = dset_session.num(station=station2)
    n = [n1, n2]

    dof.append(dset_session.meta["statistics"]["degrees of freedom"])
    variance_factor.append(dset_session.meta["statistics"]["variance factor"])
    
    if args.plot_trop:
        # Zenith wet delay
        trop_wet1 = get_state_values(dset_session, f"trop_wet-{station1}")
        trop_wet2 = get_state_values(dset_session, f"trop_wet-{station2}")
        trop_wet = np.stack((trop_wet1, trop_wet2))

        # East troposphere gradient
        trop_grad_east1 = get_state_values(dset_session, f"trop_grad-{station1}_east")
        trop_grad_east2 = get_state_values(dset_session, f"trop_grad-{station2}_east")
        trop_grad_east = np.stack((trop_grad_east1, trop_grad_east2))

        # North troposphere gradient
        trop_grad_north1 = get_state_values(dset_session, f"trop_grad-{station1}_north")
        trop_grad_north2 = get_state_values(dset_session, f"trop_grad-{station2}_north")
        trop_grad_north = np.stack((trop_grad_north1, trop_grad_north2))

        plot_sta_param(dset_session.time.utc.datetime, trop_wet, n, [f"{station1}", f"{station2}"], f"Zenith wet delay {rundate.date()} {session_code}", " Zenith Wet Delay [m]")
        plot_sta_param(dset_session.time.utc.datetime, trop_grad_east, n, [f"{station1}", f"{station2}"], f"Gradient East {rundate.date()} {session_code}", "Gradient East [m]")
        plot_sta_param(dset_session.time.utc.datetime, trop_grad_north, n, [f"{station1}", f"{station2}"], f"Gradient North {rundate.date()} {session_code}", "Gradient North [m]")


    # Collect baseline information from each session
    bl = dset_session.meta.get(baseline1) or dset_session.meta.get(baseline2)

    if bl:
        baseline_dates.append(rundate)
        baseline_length.append(bl['baseline_length'])
        baseline_length_ferr.append(bl['baseline_length_ferr'])
        num_obs_bs.append(dset_session.num(baseline=baseline1) + dset_session.num(baseline=baseline2))

    # Collect rms information for each session
    rms.append(dset_session.rms("residual"))
    rms_sta_1.append(dset_session.rms("residual", station=station1))
    rms_sta_2.append(dset_session.rms("residual", station=station2))

    if args.plot_residuals:
        plot_residuals(dset_session, station1, station2)

    if args.plot_skycoverage:
        plot_skycoverage(dset_session, station1)
        plot_skycoverage(dset_session, station2)
    
# Global plots

# Get local tie info
data = np.genfromtxt("local_ties.txt", dtype=None, encoding="utf-8", names="baseline, length")
data = dict(zip(data["baseline"], data["length"]))
local_tie = data.get(baseline1) or data.get(baseline2)

plot_statistics(dates, dof, variance_factor, colors[0])
plot_baseline(baseline_dates, baseline_length, baseline_length_ferr, num_obs_bs, station1, station2, local_tie)

if args.plot_num_obs:
    plot_session_hist(num_obs_sta1, num_obs_sta2)

plot_residual_rms(dates, rms, rms_sta_1, rms_sta_2, station1, station2, colors)

# Plot timeseries of estimated station coordinate corrections
plot_station_pos(dset_ts, station1)
plot_station_pos(dset_ts, station2)
