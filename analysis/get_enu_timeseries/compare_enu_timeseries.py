import argparse

import numpy as np 
import matplotlib.pyplot as plt



parser = argparse.ArgumentParser(epilog="Example: python compare_enu_timeseries.py  NYALES20 2023a 2025a")
parser.add_argument("station",
                    help="Compare timeseries for selected station", type=str)
parser.add_argument("id", help="Id for the timeseries to compare.", nargs="+", type=str)

parser.add_argument("--save_fig", help="Save plot as PNG", action="store_true")
args = parser.parse_args()

station = args.station
ids = args.id
data = {}

fig, axs = plt.subplots(3)
plt.suptitle(station)
axs[0].set_ylabel("East [m]")
axs[1].set_ylabel("North [m]")
axs[2].set_ylabel("Up [m]")


for  _id in ids:
    filename = f"{station}_enu_estimates_{_id}.txt"
    data[_id] = np.genfromtxt(filename, names="date, decimalyear, e, n, u, fe, fn, fu, ax, ay, az")
    axs[0].errorbar(data[_id]["decimalyear"], data[_id]["e"], marker=".", yerr=data[_id]["fe"], label=_id)
    axs[1].errorbar(data[_id]["decimalyear"], data[_id]["n"], marker=".", yerr=data[_id]["fn"], label=_id)
    axs[2].errorbar(data[_id]["decimalyear"], data[_id]["u"], marker=".", yerr=data[_id]["fu"], label=_id)

plt.legend()
plt.tight_layout()

if args.save_fig:
    plt.savefig(f"{station}_enu_{'_'.join(ids)}", dpi=300)

plt.show()
