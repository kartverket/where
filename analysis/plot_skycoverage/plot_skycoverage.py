import argparse
from datetime import date, datetime
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from where.data import dataset3 as dataset
from where.lib import config
from where import apriori


def get_residuals(dset_session, station):
    idx1 = dset_session.filter(station_1=station)
    idx2 = dset_session.filter(station_2=station)
    elevation1 = dset_session.site_pos_1.elevation[idx1]
    elevation2 = dset_session.site_pos_2.elevation[idx2]
    azimuth1 = dset_session.site_pos_1.azimuth[idx1]
    azimuth2 = dset_session.site_pos_2.azimuth[idx2]
    residual1 = dset_session.residual[idx1]
    residual2 = dset_session.residual[idx2]
    residuals = np.concatenate((residual1, residual2))
    elevation = np.concatenate((elevation1, elevation2))
    azimuth = np.concatenate((azimuth1, azimuth2))

    #el1 = np.degrees(elevation1) < 15
    #el2 = np.degrees(elevation2) < 15
    #if station == "NYALE13S" and (any(el1) or any(el2)):
    #    num_el = np.sum(el1) + np.sum(el2)
    #    print(f"Elevation below 15 degrees for {num_el} obs for {station} in session {dset_session.vars['rundate']} {dset_session.vars['session_code']}")
    return list(np.abs(residuals)), list(elevation), list(azimuth)

def plot_skycoverage(residuals, elevation, azimuth, station, label):
    title = f"{station}"

    fig = plt.figure(figsize=(8, 8), dpi=150)
    plt.set_cmap("cool")

    norm = mpl.colors.Normalize(vmin=0, vmax=0.1)
    ax = fig.add_subplot(111, projection="polar")
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    elevation = np.degrees(elevation)
    # Invert data because inverting axis is hard
    elevation = 90 - elevation
    im = ax.scatter(azimuth, elevation, c=residuals, marker=".", alpha=0.75, norm=norm)
    ax.set_title(title)
    # Change labels because elevation is "inverted"
    degree_sign = u'\N{DEGREE SIGN}'
    r_labels = [
        '90' + degree_sign,
        '',
        '60' + degree_sign,
        '',
        '30' + degree_sign,
        '',
        '0' + degree_sign,
    ]
    ax.set_rgrids(range(1, 106, 15), r_labels)
    theta_labels = [
        'N \n 0' + degree_sign, 
        '45' + degree_sign,
        'E \n 90' + degree_sign,
        '135' + degree_sign,
        'S \n 180' + degree_sign,
        '225' + degree_sign,
        'W \n 270' + degree_sign,
        '315' + degree_sign,
    ]
    ax.set_thetagrids(range(0, 360, 45), theta_labels)
    cbar = plt.colorbar(im, use_gridspec=True, pad=0.1)
    cbar.set_label("Absolute value of residual after estimation [m]")
    os.makedirs(f"img/{dset_id}/", exist_ok=True)
    fig.savefig(f"img/{dset_id}/Skyplot_{station}_{label}_{dset_id}.png", bbox_inches="tight")
    plt.close()



# Program starts execution here
    

parser = argparse.ArgumentParser(epilog="Example: python plot_skycoverage.py --id nyales --station NYALE13S",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--id", help="Dataset id of result files.", type=str, default="")
parser.add_argument("--start", help="Start date to look for sessions in master files. Format: YYYY-mm-dd", type=date.fromisoformat, default=date.min)
parser.add_argument("--end", help="End date to look for sessions in master files. Format:YYYY-mm-dd", type=date.fromisoformat, default=date.max)
parser.add_argument("--station", help="Name of a station", nargs="+", type=str, default="NYALE13S")
parser.add_argument("--session_wise", help="Enable this to get one plot per session. This is time consuming if the time window is large", action="store_true", default=False)
args = parser.parse_args()

dset_id = args.id
stations = list(args.station) if isinstance(args.station, str) else args.station
start = args.start
end = args.end

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
dates = np.array([datetime.strptime(dt, "%Y-%m-%d").date() for dt in dset_ts.rundate[idx]])
idx2 = np.logical_and(dates >= start, dates <= end)
session_codes = dset_ts.session_code[idx][idx2]
num_obs_total = np.sum(dset_ts.num_obs_estimate[idx])
print(f"{dset_id}: Number of observations total: {num_obs_total}")
# Accumulated residuals, elevation and azimuth angles
acc_r = {}
acc_e = {}
acc_a = {}
for s in stations:
    acc_r[s] = []
    acc_e[s] = []
    acc_a[s] = []

# Plot per session
for i, (rundate, session_code) in enumerate(zip(dates[idx2], session_codes)):
    print(f"{rundate} {session_code}")

    dset_session = dataset.Dataset.read(
        rundate=rundate, pipeline=pipeline, stage="postprocess", label="0", session_code=session_code, id=dset_id
    )
    for s in stations:
        r, e, a = get_residuals(dset_session, s)
        acc_r[s] += r
        acc_e[s] += e
        acc_a[s] += a
        if args.session_wise and r:
            plot_skycoverage(r, e, a, s, label=f"{rundate}_{session_code}")

for s in stations:
    plot_skycoverage(acc_r[s], acc_e[s], acc_a[s], s, label=f"{start}_{end}")
    print(f"{dset_id}: Number of observations total for {s}: {len(acc_r[s])}")

