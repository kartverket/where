import argparse
import calendar
from datetime import date, datetime, timedelta
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
import matplotlib.dates as mdates


from where.data import dataset3 as dataset
from where.lib import config


# Program starts execution here
parser = argparse.ArgumentParser(epilog="Example: python plot_num_obs.py --id nyale13s0 --stations NYALE13S NYALES20 NYALE13N",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--id", help="Dataset id of result files.", type=str, default="nyale13s0")
parser.add_argument("--stations", help="Name of the stations", nargs="+", type=str, default=["NYALE13S", "NYALES20", "NYALE13N"])
parser.add_argument("--normalized", help="Enable this flag to plot normalized bar plots.", action="store_true")
parser.add_argument("--combined", help="Enable this flag to plot all stations in one plots.", action="store_true")
parser.add_argument("--export_to_csv", help="Save data to one csv file per station.", action="store_true")
parser.add_argument("--years", help="Which year to create plots for. Default is all years found in result files.", nargs="+", type=int)
parser.add_argument("--same_scale", help="Enable this flag to force all plots to use the same scale.", action="store_true")
parser.add_argument("--stack", help="Enable this flag to stack all months on from different years", action="store_true")
parser.add_argument("--extension", help="File extension for saved images", default="{ext}")
args = parser.parse_args()

dset_id = args.id
ext = args.extension
stations = sorted(args.stations)
pipeline = "vlbi"

img_dir = f"img/{dset_id}"
csv_dir = f"csv/{dset_id}"
os.makedirs(img_dir, exist_ok=True)
os.makedirs(csv_dir, exist_ok=True)

# Setup config to make apriori modules work
config.set_analysis(rundate=None, pipeline=pipeline)
config.files.profiles = [pipeline]
config.read_pipeline(pipeline)

# Get data from timeseries dataset
dset_ts = dataset.Dataset.read(
    rundate=date(1970, 1, 1), pipeline=pipeline, stage="timeseries", label="0", session_code="", id=dset_id
)

num_sta = len(stations)
min_year = 9999
max_year = 0
max_ylim = 0
data = {}

for station in stations:
    sta_data = {}
    idx = dset_ts.filter(station=station)
    sta_data["date"] = np.array([datetime.strptime(dt, "%Y-%m-%d") for dt in dset_ts.rundate[idx]])
    sta_data["session_code"] = dset_ts.session_code[idx]
    sta_data["session_code"] = dset_ts.session_code[idx]
    sta_data["scheduled"] = dset_ts.num_obs_schedule[idx]
    sta_data["correlated"] = dset_ts.num_obs_read[idx]
    sta_data["used"] = dset_ts.num_obs_estimate[idx]
    data[station] = sta_data
    
    max_scheduled = np.nan_to_num(sta_data["scheduled"]).max()
    max_ylim = max_scheduled if max_scheduled > max_ylim else max_ylim
    min_year_sta = sta_data["date"].min().year
    max_year_sta = sta_data["date"].max().year
    min_year = min_year_sta if min_year_sta < min_year else min_year
    max_year = max_year_sta if max_year_sta > max_year else max_year

ylim = (0, max_ylim + max_ylim * 0.01)

if args.years:
    years = args.years
else:
    years = list(range(min_year, max_year+1))

# Used if args.stack is set    
stacked = {}
MONTHS = list(calendar.month_abbr)[1:]

for year in years:
    
    if args.combined:
        # One plot per year with all stations in the same plot, values are absolute
        fig_height = 4 * num_sta
        fig, axs = plt.subplots(num_sta, figsize=(12, fig_height), dpi=150, sharex=True)
        for i, station in enumerate(stations):
            idx = np.logical_and(data[station]["date"] <= datetime(year, 12, 31), data[station]["date"] >= datetime(year, 1, 1))
            plot_data = data[station]
            axs[i].bar(plot_data["date"][idx], plot_data["scheduled"][idx], label="Scheduled", width=2)
            axs[i].bar(plot_data["date"][idx], plot_data["correlated"][idx], label="Correlated", width=2)
            axs[i].bar(plot_data["date"][idx], plot_data["used"][idx], label="Used", width=2)
            xlim = (datetime(year, 1, 1), datetime(year, 12, 31))
            axs[i].set_xlim(xlim)
            if args.same_scale:
                axs[i].set_ylim(ylim)
            axs[i].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
            axs[i].xaxis.set_major_locator(mt.LinearLocator(numticks=7))
            axs[i].set_ylabel(station)
            axs[i].legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
            fig.autofmt_xdate()
    
        fig.suptitle(f"Absolute number of observations for {year}")
        fig.subplots_adjust(top=0.95)
        plt.savefig(f"{img_dir}/Num_obs_{'_'.join(stations)}_{year}_{dset_id}.{ext}", bbox_inches="tight")
        plt.close()

    if args.combined and args.normalized:
        # One plot per year with all stations in same plot, values are normalized
        fig_height = 4 * num_sta
        fig, axs = plt.subplots(num_sta, figsize=(12, fig_height), dpi=150, sharex=True)
        for i, station in enumerate(stations):
            idx = np.logical_and(data[station]["date"] <= datetime(year, 12, 31), data[station]["date"] >= datetime(year, 1, 1))
            plot_data = data[station]
            norm_factor = 100/plot_data["scheduled"][idx]
            axs[i].bar(plot_data["date"][idx], plot_data["scheduled"][idx]*norm_factor, label="Scheduled", width=2)
            axs[i].bar(plot_data["date"][idx], plot_data["correlated"][idx]*norm_factor, label="Correlated", width=2)
            axs[i].bar(plot_data["date"][idx], plot_data["used"][idx]*norm_factor, label="Used", width=2)
            xlim = (datetime(year, 1, 1), datetime(year, 12, 31))
            axs[i].set_xlim(xlim)
            #if args.same_scale:
            #    axs[i].set_ylim(ylim)

            axs[i].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
            axs[i].xaxis.set_major_locator(mt.LinearLocator(numticks=7))
            axs[i].set_ylabel(station)
            axs[i].legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
            fig.autofmt_xdate()
    
        fig.suptitle(f"Normalized number of observations for {year}")
        fig.subplots_adjust(top=0.95)
        plt.savefig(f"{img_dir}/Num_obs_{'_'.join(stations)}_normalized_{year}_{dset_id}.{ext}", bbox_inches="tight")
        plt.close()

    # Default. Always plot at least these plots
    # One plot per year and one plot per station, values are absolute
    for i, station in enumerate(stations):
        idx = np.logical_and(data[station]["date"] <= datetime(year, 12, 31), data[station]["date"] >= datetime(year, 1, 1))
        if not any(idx):
            continue
        fig = plt.figure(figsize=(12,8), dpi=150)
        plot_data = data[station]
        plt.bar(plot_data["date"][idx], plot_data["scheduled"][idx], label="Scheduled", width=2)
        plt.bar(plot_data["date"][idx], plot_data["correlated"][idx], label="Correlated", width=2)
        plt.bar(plot_data["date"][idx], plot_data["used"][idx], label="Used", width=2)
        xlim = (datetime(year, 1, 1), datetime(year, 12, 31))
        plt.xlim(xlim)
        if args.same_scale:
            plt.ylim(ylim)

        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        plt.gca().xaxis.set_major_locator(mt.LinearLocator(numticks=7))
        plt.ylabel(station)
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        fig.autofmt_xdate()

        fig.suptitle(f"Absolute number of observations for {year}")
        #fig.subplots_adjust(top=0.95)
        plt.savefig(f"{img_dir}/Num_obs_{station}_{year}_{dset_id}.{ext}", bbox_inches="tight")
        plt.close()

    if args.normalized:
        # One plot per year and one plot per station, values are normalized
        for i, station in enumerate(stations):
            idx = np.logical_and(data[station]["date"] <= datetime(year, 12, 31), data[station]["date"] >= datetime(year, 1, 1))
            if not any(idx):
                continue
            fig = plt.figure(figsize=(12,8), dpi=150)
            plot_data = data[station]
            norm_factor = 100/plot_data["scheduled"][idx]
            plt.bar(plot_data["date"][idx], plot_data["scheduled"][idx]*norm_factor, label="Scheduled", width=2)
            plt.bar(plot_data["date"][idx], plot_data["correlated"][idx]*norm_factor, label="Correlated", width=2)
            plt.bar(plot_data["date"][idx], plot_data["used"][idx]*norm_factor, label="Used", width=2)
            xlim = (datetime(year, 1, 1), datetime(year, 12, 31))
            plt.xlim(xlim)
            #if args.same_scale:
            #    plt.ylim(ylim)

            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
            plt.gca().xaxis.set_major_locator(mt.LinearLocator(numticks=7))
            plt.ylabel(station)
            plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
            fig.autofmt_xdate()
    
            plt.suptitle(f"Normalized number of observations for {year}")
            #fig.subplots_adjust(top=0.95)
            plt.savefig(f"{img_dir}/Num_obs_{station}_normalized_{year}_{dset_id}.{ext}", bbox_inches="tight")
            plt.close()

    if args.stack:
        for station in stations:
            sta_dict = stacked.setdefault(station, {})
            sta_dict.setdefault("scheduled", np.zeros(len(MONTHS)))
            sta_dict.setdefault("correlated", np.zeros(len(MONTHS)))
            sta_dict.setdefault("used", np.zeros(len(MONTHS)))
            for i, m in enumerate(MONTHS):
                #month_dict = sta_dict.setdefault(m, {})
                start = datetime(year, i+1, 1)
                end = datetime(year, i+1, calendar.monthrange(year, i+1)[1])
                idx = np.logical_and(data[station]["date"] >= start, data[station]["date"] <= end)
                #month_dict.setdefault("scheduled", 0)
                #month_dict.setdefault("correlated", 0)
                #month_dict.setdefault("used", 0)
                #import IPython; IPython.embed()
                # Ignore sessions without information about shedule
                idx_nan = np.isnan(data[station]["scheduled"][idx])
                sta_dict["scheduled"][i] += np.nansum(data[station]["scheduled"][idx][~idx_nan])
                sta_dict["correlated"][i] += np.nansum(data[station]["correlated"][idx][~idx_nan])
                sta_dict["used"][i] += np.nansum(data[station]["used"][idx][~idx_nan])
    
    
        
    if args.export_to_csv:
        for station in stations:
            with open(f"{csv_dir}/num_obs_{station}.csv", "w") as fid:
                fid.write(f"session code, date, scheduled, correlated, used\n")
                scheduled = np.nan_to_num(data[station]["scheduled"])
                correlated = np.nan_to_num(data[station]["correlated"])
                used = np.nan_to_num(data[station]["used"])
                session_codes = data[station]["session_code"]
                dates = data[station]["date"]
                for sc, d, s, r, u in zip(session_codes, dates, scheduled, correlated, used): 
                    fid.write(f"{sc:8}, {d:%Y-%m-%d}, {int(s)}, {int(r)}, {int(u)} \n")

if args.stack:
    for station in stacked.keys():
        #import IPython; IPython.embed()
        plt.bar(MONTHS, stacked[station]["scheduled"], label="Scheduled")
        plt.bar(MONTHS, stacked[station]["correlated"], label="Correlated")
        plt.bar(MONTHS, stacked[station]["used"], label="Used")
        plt.legend()
        plt.ylabel("Number of observations per month")
        plt.title(f"{station} from {min_year} to {max_year}")
        plt.savefig(f"{img_dir}/Num_obs_per_month_{station}_{min_year}_{max_year}_{dset_id}.{ext}", bbox_inches="tight")
        plt.close()
        
        # Normalized 
        norm_factor = 100/stacked[station]["scheduled"]
        plt.bar(MONTHS, stacked[station]["scheduled"]*norm_factor, label="Scheduled")
        plt.bar(MONTHS, stacked[station]["correlated"]*norm_factor, label="Correlated")
        plt.bar(MONTHS, stacked[station]["used"]*norm_factor, label="Used")
        plt.legend()
        plt.ylabel("Normalized number of observations per month")
        plt.title(f"{station} from {min_year} to {max_year}")
        plt.savefig(f"{img_dir}/Num_obs_per_month_{station}_normalized_{min_year}_{max_year}_{dset_id}.{ext}", bbox_inches="tight")
        plt.close()
        
        