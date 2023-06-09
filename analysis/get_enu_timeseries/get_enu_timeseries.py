import argparse
from datetime import datetime, date

import numpy as np
import matplotlib.pyplot as plt

from where import apriori
from where.lib import config
from where.lib import rotation
from where.data import dataset3 as dataset
from where.data.time import Time
from where.data.position import PositionDelta

pipeline = "vlbi"

parser = argparse.ArgumentParser(epilog="Example: python get_enu_timeseries.py --station NYALES20 --id nnstest")
parser.add_argument("--station", 
                    help="Create timeseries for selected station(s). If omitted, timeseries for all stations will be created",
                    nargs="*",
                    type=str)
parser.add_argument("--id", help="Dataset id for timeseries dataset.", type=str, default="")
parser.add_argument("--show_fig", help="Show plot of ENU-timeseries", action="store_true")

args = parser.parse_args()
stations = args.station
show_fig = args.show_fig
dset_id = args.id


# Setup config to make apriori modules work
config.set_analysis(rundate=None, pipeline=pipeline)
config.files.profiles = [pipeline]
config.read_pipeline(pipeline)

# Get data from timeseries dataset
dset_ts = dataset.Dataset.read(rundate=date(1970,1,1), pipeline=pipeline, stage="timeseries", label="0", session_code="", id=dset_id)
if stations is None:
    stations = list(dset_ts.unique("station"))
    stations.remove("all")

for station in stations:
    idx = dset_ts.filter(station=station, status="unchecked")
    t = dset_ts.time[idx]
    if len(t) <= 1:
        print(f"Skipping {station} due to only 1 data point")
        continue
    # Sort data based on time
    sort_idx = np.argsort(t.mjd)
    t = t[sort_idx]
    trf = apriori.get("trf", time=t, reference_frames="itrf:2020")
    names = apriori.get("vlbi_station_codes")

    if names[station]["cdp"] not in trf:
        print(f"Skipping {station} since it is not in {', '.join(trf.reference_frames)}")
        continue

    pos = trf[names[station]["cdp"]].pos
    lat, lon, _ = pos.llh.T
    trs2enu = rotation.enu2trs(lat, lon)
    enu2trs = rotation.trs2enu(lat, lon)

    dpos = PositionDelta(val=np.squeeze(dset_ts.neq_vlbi_site_pos[idx][sort_idx]), system="trs", ref_pos=pos)
    dpos_cov_xyz = dset_ts.neq_vlbi_site_pos_cov_[idx][sort_idx]
    dpos_ferr_xyz = np.sqrt(dpos_cov_xyz.diagonal(axis1=1, axis2=2))
    dpos_cov_enu = trs2enu @ dpos_cov_xyz @ enu2trs
    dpos_ferr_enu = np.sqrt(dpos_cov_enu.diagonal(axis1=1, axis2=2))
    
    #print(f"Max ferr east : {t[np.argmax(dpos_ferr_enu.T[0])].datetime}, {dpos_ferr_enu.T[0][np.argmax(dpos_ferr_enu.T[0])]:8.4f}")
    #print(f"Max ferr north: {t[np.argmax(dpos_ferr_enu.T[1])].datetime}, {dpos_ferr_enu.T[1][np.argmax(dpos_ferr_enu.T[1])]:8.4f}")
    #print(f"Max ferr up   : {t[np.argmax(dpos_ferr_enu.T[2])].datetime}, {dpos_ferr_enu.T[2][np.argmax(dpos_ferr_enu.T[2])]:8.4f}")
    #print(f"Max east : {t[np.argmax(dpos.enu.east)].datetime}, {dpos.x[np.argmax(dpos.enu.east)]:8.4f}")
    #print(f"Max north: {t[np.argmax(dpos.enu.north)].datetime}, {dpos.x[np.argmax(dpos.enu.north)]:8.4f}")
    #print(f"Max up   : {t[np.argmax(dpos.enu.up)].datetime}, {dpos.x[np.argmax(dpos.enu.up)]:8.4f}")

    header = f"""# East North and Up estimates for {station} produced by Where
# Format:
#
# date :  date for decimal year epoch
# decimal year: given at the mean epoch of session
# estimated correction in east: [m]
# estimated correction in north: [m]
# estimated correction in up: [m]
# formal error of east estimate: [m]
# formal error of north estimate: [m]
# formal error of up estimate: [m]
# a priori position in x coordinate: [m] (ITRF2020)
# a priori position in y coordinate: [m] (ITRF2020)
# a priori position in z coordinate: [m] (ITRF2020)
"""

    outfile=f"{station}_enu_estimates_{dset_id}.txt"
    #print(f"Writing {outfile}")
    with open(outfile, "w") as fid:
        fid.write(header)
        for epoch, apriori_pos, delta_pos, ferr in zip(t, pos, dpos.enu, dpos_ferr_enu):
            fid.write(f"{epoch.date} {epoch.decimalyear:17.12f} {delta_pos.east: 7.4f} {delta_pos.north: 7.4f} {delta_pos.up: 7.4f} {ferr[0]:7.4f} {ferr[1]:7.4f} {ferr[2]:7.4f} {apriori_pos.x: 12.4f} {apriori_pos.y: 12.4f} {apriori_pos.z: 12.4f}\n") 

    if show_fig:
        fig, axs = plt.subplots(3)
        plt.suptitle(station)
        axs[0].errorbar(t.datetime, dpos.enu.east, yerr=dpos_ferr_enu.T[0], fmt=".")
        axs[0].set_ylabel("East [m]")
        axs[1].errorbar(t.datetime, dpos.enu.north, yerr=dpos_ferr_enu.T[1], fmt=".")
        axs[1].set_ylabel("North [m]")
        axs[2].errorbar(t.datetime, dpos.enu.up, yerr=dpos_ferr_enu.T[2], fmt=".")
        axs[2].set_ylabel("Up [m]")
        plt.show()
