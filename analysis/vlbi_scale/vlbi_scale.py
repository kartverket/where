import argparse
from datetime import datetime, date, timedelta
import os

import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt

from where.lib import config
from where.data import dataset3 as dataset
from where.data.time import Time

# Setup input argument parser for script 
parser = argparse.ArgumentParser(epilog="Example: python vlbi_scale.py --year=2013.75 --id=test --rms_limit=5 --volume_limit=10 --save_fig --show_fig",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--id", help="Dataset id for the dataset timeseries", type=str, default="")
parser.add_argument("--year", help="Decimal year of discontinuity in scale", type=float, default=2013.75)
parser.add_argument("--rms_limit", help="Outliers larger than rms_limit*RMS are discarded", type=float, default=5)
parser.add_argument("--volume_limit", help="Networks smaller than the given limit (in Mm^3) are discarded", type=float, default=10)
parser.add_argument("--save_fig", help="Plot is saved as a png file", action="store_true", default=False)
parser.add_argument("--show_fig", help="Plot is showed interactively", action="store_true", default=False)
parser.add_argument("--make_csv", help="Data is saved as a csv file", action="store_true", default=False)

# Parse input args
args = parser.parse_args()
dset_id = args.id
solution_id = dset_id if dset_id else "default"
discont_epoch = Time(args.year, scale="utc", fmt="decimalyear")
rms_limit = args.rms_limit
volume_limit = args.volume_limit
save_fig = args.save_fig
show_fig = args.show_fig
make_csv = args.make_csv

# Create output directories
imgdir = "img"
csvdir = "csv"
os.makedirs(imgdir, exist_ok=True)
os.makedirs(csvdir, exist_ok=True)

# Setup config to make apriori modules work
pipeline = "vlbi"
config.set_analysis(rundate=None, pipeline=pipeline)
config.files.profiles = [pipeline]
config.read_pipeline(pipeline)

# Get data from timeseries dataset
dset_ts = dataset.Dataset.read(rundate=date(1970,1,1), pipeline=pipeline, stage="timeseries", label="0", session_code="", id=dset_id)

# Select and discard sessions
idx_volume = dset_ts.network_volume > volume_limit
idx_all_good = dset_ts.filter(station="all", status="unchecked")
idx_nan = np.isnan(dset_ts.neq_helmert_scale)
idx_1 = np.logical_and(np.logical_and(idx_volume, idx_all_good), ~idx_nan)
scale_rms = np.sqrt(np.mean(dset_ts.neq_helmert_scale[idx_1]**2))
idx_rms = np.abs(dset_ts.neq_helmert_scale) < rms_limit*scale_rms
idx_2 = np.logical_and(idx_1, idx_rms)
sort_idx_2 = np.argsort(dset_ts.rundate[idx_2])


print(f"Total number of sessions: {np.sum(idx_all_good)}")
print(f"Discaring sessions with network volume smaller than 10Mm^3")
print(f"Discaring sessions with scale larger than {rms_limit*scale_rms:6.4f} ppb.")
print(f"Discaring sessions due to computation failure")
print(f"Total number of sessions left: {np.sum(idx_2)}")

# Split data by input year
idx_before = dset_ts.time.utc.datetime <= discont_epoch.utc.datetime
idx_after = dset_ts.time.utc.datetime > discont_epoch.utc.datetime
idx_3 = np.logical_and(idx_2, idx_before)
idx_4 = np.logical_and(idx_2, idx_after)

# Perform linear regression
result_before = linregress(dset_ts.time.utc.decimalyear[idx_3], dset_ts.neq_helmert_scale[idx_3])
fitted_before = result_before.intercept + result_before.slope*dset_ts.time.utc.decimalyear[idx_3]
result_after = linregress(dset_ts.time.utc.decimalyear[idx_4], dset_ts.neq_helmert_scale[idx_4])
fitted_after = result_after.intercept + result_after.slope*dset_ts.time.utc.decimalyear[idx_4]
mean_scale_before = np.mean(dset_ts.neq_helmert_scale[idx_3])
mean_scale_after = np.mean(dset_ts.neq_helmert_scale[idx_4])
mean_time_before = dset_ts.time[idx_3].mean.utc.datetime
mean_time_after = dset_ts.time[idx_4].mean.utc.datetime
print(f"Before {discont_epoch.utc.decimalyear}: Slope: {result_before.slope: 6.4f} [ppb/year], Mean = {mean_scale_before: 6.4f} [ppb]")
print(f"After {discont_epoch.utc.decimalyear}:  Slope: {result_after.slope: 6.4f} [ppb/year], Mean = {mean_scale_after: 6.4f} [ppb]")

# Plot scale with linear regression
plt.figure(figsize=(12,8), dpi=150)
plt.scatter(dset_ts.time.utc.datetime[idx_2][sort_idx_2], dset_ts.neq_helmert_scale[idx_2][sort_idx_2], c=dset_ts.network_volume[idx_2][sort_idx_2], marker=".")
plt.axvline(discont_epoch.utc.datetime, color="black")
plt.plot(dset_ts.time.utc.datetime[idx_3], fitted_before, color="red")
plt.plot(dset_ts.time.utc.datetime[idx_4], fitted_after, color="red")
plt.gca().annotate(f"{result_before.slope: 6.4f} [ppb/year]", xy=(mean_time_before, mean_scale_before), xycoords="data",
    xytext=(0.3, 0.05), textcoords='axes fraction', arrowprops=dict(facecolor="black", width=1, headwidth=4))
plt.gca().annotate(f"{result_after.slope: 6.4f} [ppb/year]", xy=(mean_time_after, mean_scale_after), xycoords="data",
    xytext=(0.7, 0.95), textcoords='axes fraction', arrowprops=dict(facecolor="black", width=1, headwidth=4))
cbar = plt.colorbar()
cbar.set_label(f"Network volume [Mm^3]")
plt.title(f"VLBI scale with linear regression before and after {discont_epoch}")
plt.ylabel(f"Scale [ppb]")
if save_fig:
    plt.savefig(f"{imgdir}/VLBI_scale_linear_regression_{discont_epoch}_{rms_limit}_{volume_limit}_{solution_id}.png", bbox_inches="tight")
if show_fig:
    plt.show()
plt.close()


# Plot the other Helmert parameters also
fields = ["neq_helmert_T_X", "neq_helmert_T_Y", "neq_helmert_T_Z", "neq_helmert_alpha", "neq_helmert_beta", "neq_helmert_gamma"]
units = ["mm", "mm", "mm", "mas", "mas", "mas"]
for f, u in zip(fields,units):
    name = f.replace("neq_helmert_", "")
    plt.figure(figsize=(12,8), dpi=150)
    plt.scatter(dset_ts.time[idx_2][sort_idx_2].utc.datetime, dset_ts[f][idx_2][sort_idx_2], c=dset_ts.network_volume[idx_2][sort_idx_2], marker=".")
    cbar = plt.colorbar()
    cbar.set_label(f"Network volume [Mm^3]")
    plt.title(f"{name}")
    plt.ylabel(f"{name} [{u}]")
    if save_fig:
        plt.savefig(f"{imgdir}/VLBI_{name}_{rms_limit}_{volume_limit}_{solution_id}.png", bbox_inches="tight")
    if show_fig:
        plt.show()
    plt.close()


# Compute running median
days = 90
delta = timedelta(days=days/2)
start_epoch = dset_ts.time[idx_2].min.utc.datetime
end_epoch = dset_ts.time[idx_2].max.utc.datetime
running_median = np.full(len(dset_ts.time[idx_2]), fill_value=np.nan)


for i, t in enumerate(dset_ts.time[idx_2][sort_idx_2].utc.datetime):
    time_lower_limit = t - delta
    if time_lower_limit < start_epoch:
        continue

    time_upper_limit = t + delta
    if time_upper_limit > end_epoch:
        continue

    idx_time_lower = dset_ts.time[idx_2][sort_idx_2].utc.datetime > time_lower_limit
    idx_time_upper = dset_ts.time[idx_2][sort_idx_2].utc.datetime < time_upper_limit
    idx_interval = np.logical_and(idx_time_lower, idx_time_upper)
    #print(f"{np.sum(idx_interval)} sessions in interval {time_lower_limit:%Y-%m-%d}-{time_upper_limit:%Y-%m-%d}")
    running_median[i] = np.median(dset_ts.neq_helmert_scale[idx_2][sort_idx_2][idx_interval])

plt.figure(figsize=(12,8), dpi=150)
plt.plot(dset_ts.time[idx_2][sort_idx_2].utc.datetime, running_median, color="red")
plt.scatter(dset_ts.time[idx_2][sort_idx_2].utc.datetime, dset_ts.neq_helmert_scale[idx_2][sort_idx_2], c=dset_ts.network_volume[idx_2][sort_idx_2], marker=".")
cbar = plt.colorbar()
cbar.set_label(f"Network volume [Mm^3]")
plt.title(f"VLBI Scale with running median ({days} days)")
plt.ylabel(f"Scale [ppb]")
if save_fig:
    plt.savefig(f"{imgdir}/VLBI_scale_running_median_{days}_{rms_limit}_{volume_limit}_{solution_id}.png", bbox_inches="tight")
if show_fig:
    plt.show()
plt.close()    

# Make CSV file
if make_csv:
    filename = f"{csvdir}/VLBI_scale_{rms_limit}_{volume_limit}_{solution_id}.csv"
    with open(filename, "w") as fid:
        header = f"date, session code, scale [ppb], t_x [m], t_y [m], t_z [m], r_1 [mas], r_2 [mas], r_3 [mas], volume [Mm^3], stations "
        fid.write(f"{header}\n")
        rundate = dset_ts.rundate[idx_2][sort_idx_2]
        session_code = dset_ts.session_code[idx_2][sort_idx_2]
        scale = dset_ts.neq_helmert_scale[idx_2][sort_idx_2]
        t_x = dset_ts.neq_helmert_T_X[idx_2][sort_idx_2]
        t_y = dset_ts.neq_helmert_T_Y[idx_2][sort_idx_2]
        t_z = dset_ts.neq_helmert_T_Z[idx_2][sort_idx_2]
        r_1 = dset_ts.neq_helmert_alpha[idx_2][sort_idx_2]
        r_2 = dset_ts.neq_helmert_beta[idx_2][sort_idx_2]
        r_3 = dset_ts.neq_helmert_gamma[idx_2][sort_idx_2]
        volume = dset_ts.network_volume[idx_2][sort_idx_2]
        for rd, sc, s, tx, ty, tz, r1, r2, r3, v in zip(rundate, session_code, scale, t_x, t_y, t_z, r_1, r_2, r_3, volume):
            session_idx = dset_ts.filter(rundate=rd)
            station_idx = dset_ts.num_obs_estimate[session_idx] > 0
            stations = "/".join(list(dset_ts.station[session_idx][station_idx])[1:]) # First station entry is 'all' for each rundate
            fid.write(f"{rd}, {sc}, {s}, {tx}, {ty}, {tz}, {r1}, {r2}, {r3}, {v}, {stations}\n")
