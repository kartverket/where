import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mt
import numpy as np
import sys
from datetime import datetime

try:
    dset_id = sys.argv[1]
except IndexError:
    dset_id = "nyale13s0"

stations = ["NYALES20", "NYALE13S", "NYALE13N"]
num_sta = len(stations)
data = {}
min_year = 9999
max_year = 0
str2date = lambda x : datetime.strptime(x.decode("utf-8"), "%Y-%m-%d")
for station in stations:
    data[station] = np.genfromtxt(f"session_stats_{station}.txt", names=True, autostrip=True, delimiter=",", dtype="U8, O, i8, i8, i8", converters={1:str2date})
    min_year_sta = data[station]["date"].min().year
    max_year_sta = data[station]["date"].max().year
    min_year = min_year_sta if min_year_sta < min_year else min_year
    max_year = max_year_sta if max_year_sta > max_year else max_year

years = list(range(min_year, max_year+1))

for year in years:
    # One plot per year with all stations in the same plot, values are absolute
    fig_height = 4 * num_sta
    fig, axs = plt.subplots(num_sta, figsize=(12, fig_height), dpi=150, sharex=True)
    for i, station in enumerate(stations):
        idx = np.logical_and(data[station]["date"] <= datetime(year, 12, 31), data[station]["date"] >= datetime(year, 1, 1))
        plot_data = data[station][idx]
        axs[i].bar(plot_data["date"], plot_data["scheduled"], label="Scheduled", width=3)
        axs[i].bar(plot_data["date"], plot_data["recovered"], label="Correlated", width=3)
        axs[i].bar(plot_data["date"], plot_data["used"], label="Used", width=3)
        xlim = (datetime(year, 1, 1), datetime(year, 12, 31))
        axs[i].set_xlim(xlim)
        axs[i].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        axs[i].xaxis.set_major_locator(mt.LinearLocator(numticks=7))
        axs[i].set_ylabel(station)
        axs[i].legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        fig.autofmt_xdate()

    fig.suptitle(f"Absolute number of observations for {year}")
    fig.subplots_adjust(top=0.95)
    plt.savefig(f"img/{dset_id}/Num_obs_{'_'.join(stations)}_{year}_{dset_id}.png", bbox_inches="tight")
    plt.close()

    # One plot per year with all stations in same plot, values are normalized
    fig_height = 4 * num_sta
    fig, axs = plt.subplots(num_sta, figsize=(12, fig_height), dpi=150, sharex=True)
    for i, station in enumerate(stations):
        idx = np.logical_and(data[station]["date"] <= datetime(year, 12, 31), data[station]["date"] >= datetime(year, 1, 1))
        plot_data = data[station][idx]
        norm_factor = 100/plot_data["scheduled"]
        axs[i].bar(plot_data["date"], plot_data["scheduled"]*norm_factor, label="Scheduled", width=3)
        axs[i].bar(plot_data["date"], plot_data["recovered"]*norm_factor, label="Correlated", width=3)
        axs[i].bar(plot_data["date"], plot_data["used"]*norm_factor, label="Used", width=3)
        xlim = (datetime(year, 1, 1), datetime(year, 12, 31))
        axs[i].set_xlim(xlim)
        axs[i].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        axs[i].xaxis.set_major_locator(mt.LinearLocator(numticks=7))
        axs[i].set_ylabel(station)
        axs[i].legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        fig.autofmt_xdate()

    fig.suptitle(f"Normalized number of observations for {year}")
    fig.subplots_adjust(top=0.95)
    plt.savefig(f"img/{dset_id}/Num_obs_{'_'.join(stations)}_normalized_{year}_{dset_id}.png", bbox_inches="tight")
    plt.close()

    # One plot per year and one plot per station, values are absolute
    #fig, axs = plt.subplots(2, figsize=(12, 8), dpi=150, sharex=True)
    for i, station in enumerate(stations):
        plt.figure(figsize=(12,8), dpi=150)
        idx = np.logical_and(data[station]["date"] <= datetime(year, 12, 31), data[station]["date"] >= datetime(year, 1, 1))
        plot_data = data[station][idx]
        plt.bar(plot_data["date"], plot_data["scheduled"], label="Scheduled", width=3)
        plt.bar(plot_data["date"], plot_data["recovered"], label="Correlated", width=3)
        plt.bar(plot_data["date"], plot_data["used"], label="Used", width=3)
        xlim = (datetime(year, 1, 1), datetime(year, 12, 31))
        plt.xlim(xlim)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        plt.gca().xaxis.set_major_locator(mt.LinearLocator(numticks=7))
        plt.ylabel(station)
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        fig.autofmt_xdate()

        fig.suptitle(f"Absolute number of observations for {year}")
        #fig.subplots_adjust(top=0.95)
        plt.savefig(f"img/{dset_id}/Num_obs_{station}_{year}_{dset_id}.png", bbox_inches="tight")
        plt.close()

    # One plot per year and one plot per station, values are normalized
    #fig, axs = plt.subplots(2, figsize=(12, 8), dpi=150, sharex=True)
    for i, station in enumerate(stations):
        plt.figure(figsize=(12,8), dpi=150)
        idx = np.logical_and(data[station]["date"] <= datetime(year, 12, 31), data[station]["date"] >= datetime(year, 1, 1))
        plot_data = data[station][idx]
        norm_factor = 100/plot_data["scheduled"]
        plt.bar(plot_data["date"], plot_data["scheduled"]*norm_factor, label="Scheduled", width=3)
        plt.bar(plot_data["date"], plot_data["recovered"]*norm_factor, label="Correlated", width=3)
        plt.bar(plot_data["date"], plot_data["used"]*norm_factor, label="Used", width=3)
        xlim = (datetime(year, 1, 1), datetime(year, 12, 31))
        plt.xlim(xlim)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        plt.gca().xaxis.set_major_locator(mt.LinearLocator(numticks=7))
        plt.ylabel(station)
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        fig.autofmt_xdate()

        plt.suptitle(f"Normalized number of observations for {year}")
        #fig.subplots_adjust(top=0.95)
        plt.savefig(f"img/{dset_id}/Num_obs_{station}_normalized_{year}_{dset_id}.png", bbox_inches="tight")
        plt.close()

