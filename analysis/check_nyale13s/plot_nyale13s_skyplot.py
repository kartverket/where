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

from midgard.math.unit import Unit
from midgard.math.constant import constant

from where.data import dataset3 as dataset
from where.data.position import PositionDelta
from where.data.time import Time
from where.lib import config
from where.lib import rotation
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

    return list(np.abs(residuals)), list(elevation), list(azimuth)

def plot_skycoverage(residuals, elevation, azimuth, station):
    title = f"{station}"

    fig = plt.figure(figsize=(8, 8), dpi=150)
    ax = fig.add_subplot(111, projection="polar")
    im = ax.scatter(azimuth, elevation, c=residuals, marker=".", alpha=0.75)
    ax.set_title(title)
    cbar = plt.colorbar(im, use_gridspec=True)
    cbar.set_label("Absolute value of residual after estimation [m]")
    plt.tight_layout()
    sub_dir = "Skyplot"
    os.makedirs(f"img/{dset_id}/{sub_dir}", exist_ok=True)
    fig.savefig(f"img/{dset_id}/{sub_dir}/{sub_dir}_{title}_{dset_id}.png")
    plt.close()



# Program starts execution here
    

parser = argparse.ArgumentParser(epilog="Example: python plot_nyale13s.py --id nyale13s0 --station NYALE13S")
parser.add_argument("--id", help="Dataset id of result files.", type=str, default="nyale13s0")
parser.add_argument("--station", help="Name of a station", nargs="+", type=str, default="NYALE13S")
args = parser.parse_args()

dset_id = args.id
stations = args.station

pipeline = "vlbi"

os.makedirs(f"img/{dset_id}", exist_ok=True)

# Setup config to make apriori modules work
config.set_analysis(rundate=None, pipeline=pipeline)
config.files.profiles = [pipeline]
config.read_pipeline(pipeline)

# Get data from timeseries dataset
dset_ts = dataset.Dataset.read(
    rundate=date(1970, 1, 1), pipeline=pipeline, stage="timeseries", label="0", session="", id=dset_id
)

idx = dset_ts.filter(station="all")
dates = np.array([datetime.strptime(dt, "%Y-%m-%d") for dt in dset_ts.rundate[idx]])
sessions = dset_ts.session[idx]

# Accumulated residuals, elevation and azimuth angles
acc_r = {}
acc_e = {}
acc_a = {}
for s in stations:
    acc_r[s] = []
    acc_e[s] = []
    acc_a[s] = []

# Plot per session
for i, (rundate, session) in enumerate(zip(dates, sessions)):
    print(f"{rundate} {session}")

    dset_session = dataset.Dataset.read(
        rundate=rundate, pipeline=pipeline, stage="postprocess", label="0", session=session, id=dset_id
    )
    for s in stations:
        r, e, a = get_residuals(dset_session, s)
        acc_r[s] += r
        acc_e[s] += e
        acc_a[s] += a

for s in stations:
    plot_skycoverage(acc_r[s], acc_e[s], acc_a[s], s)

