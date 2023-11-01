import argparse
from datetime import date, datetime, timedelta
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
import matplotlib.dates as mdates


from where.data import dataset3 as dataset
from where.lib import config


def make_barplot(data, station):
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


# Program starts execution here
parser = argparse.ArgumentParser(epilog="Example: python plot_num_obs.py --id nyale13s0 --stations NYALE13S NYALES20 NYALE13N",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--id", help="Dataset id of result files.", type=str, default="nyale13s0")
parser.add_argument("--stations", help="Name of the stations", nargs="+", type=str, default=["NYALE13S", "NYALES20", "NYALE13N"])
args = parser.parse_args()

dset_id = args.id
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

for station in args.stations:
    num_obs_sta = {}
    idx = dset_ts.filter(station=station)
    num_obs_sta["date"] = np.array([datetime.strptime(dt, "%Y-%m-%d") for dt in dset_ts.rundate[idx]])
    num_obs_sta["session_code"] = dset_ts.session_code[idx]
    num_obs_sta["schedule"] = dset_ts.num_obs_schedule[idx]
    num_obs_sta["read"] = dset_ts.num_obs_read[idx]
    num_obs_sta["estimate"] = dset_ts.num_obs_estimate[idx]
    make_barplot(num_obs_sta, station)

