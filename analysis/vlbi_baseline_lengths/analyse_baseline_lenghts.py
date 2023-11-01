import argparse
from datetime import datetime
import glob
import cartopy.crs as ccrs
from matplotlib import cm
from matplotlib import colors
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mt
import os

import numpy as np

from midgard.math.unit import Unit

from where.data.time import Time
from where.lib import config
from where import parsers

DATE_FMT = "%Y-%m-%d"

def read_session_baseline(dset_id=None):
    """Reads baseline files from session work directories.

    Args:
        dset_id:    Dataset id.

    Returns:
        Dictionary with data from files. Key: Baseline name. Value: Dictionary with data
    """
    print(f"Reading baseline session files from Where work directory")

    file_vars = config.create_file_vars(rundate=None, pipeline="vlbi", id=dset_id)
    dates = sorted(config.files.glob_variable("vlbi_baseline_lengths", "date", r"[0-9]{8}", file_vars=file_vars))

    baselines = {}
    for d in dates:
        rundate = datetime.strptime(d, "%Y%m%d")
        file_vars.update(date=d)
        session_codes = sorted(
            config.files.glob_variable("vlbi_baseline_lengths", "session_code", r"[_a-zA-Z0-9]*", file_vars=file_vars)
        )
        for session_code in session_codes:
            session_vars = config.create_file_vars(rundate.date(), "vlbi", session_code=session_code, id=dset_id)
            p = parsers.parse_key("vlbi_baseline_lengths", file_vars=session_vars, parser_name="vlbi_baseline_lengths")
            if p.data:
                for baseline, value in p.data.items():
                    bl = baselines.setdefault(baseline, {})
                    bl.setdefault("date", []).append(rundate)
                    bl.setdefault("session_code", []).append(session_code)
                    bl.setdefault("length", []).append(value["length"])
                    bl.setdefault("ferr", []).append(value["ferr"])
                    bl.setdefault("num_obs", []).append(value["num_obs"])
                    bl.setdefault("lat_1", []).append(value["lat_1"])
                    bl.setdefault("lat_2", []).append(value["lat_2"])
                    bl.setdefault("lon_1", []).append(value["lon_1"])
                    bl.setdefault("lon_2", []).append(value["lon_2"])
    # Convert data to np.array
    return {k1:{k2:np.array(v2) for k2, v2 in v1.items()} for k1,v1 in baselines.items()}


def read_baselines(id_txt="default"):
    """Read baseline files from 'txt' directories.

    Args:
        id_txt:    String with dataset id

    Returns:
        Dictionary with data from files. Key: Baseline name. Value: Dictionary with data
    """
    print(f"Reading baseline files in txt/{id_txt}")
    
    files = glob.glob(f"txt/{id_txt}/Baseline_*_{id_txt}.txt")
    files.remove(f"txt/{id_txt}/Baseline_length_repeatability_{id_txt}.txt")
    baselines = {}
    for bl_file in files:
        with open(bl_file) as fid:
            baseline = fid.readline().strip("#").strip()
        
        str2date = lambda x: datetime.strptime(x.decode('utf-8'), DATE_FMT)
        #print(f"Reading {bl_file}")
        data = np.genfromtxt(bl_file, names="date, session_code, num_obs, length, ferr", dtype="object, U7, i8, f8, f8", converters={0: str2date})
        data = np.atleast_1d(data)
        baselines[baseline] = {k:np.array(data[k]) for k in data.dtype.names}
    return baselines


def save_baselines(baselines, id_txt):
    """Save baseline length for each baseline in a text file.

    Saves the time evolution of the baseline length with formal errors. Each baseline is saved to
    individual files in the 'txt' directory.

    Args:
        baselines:        Dictionary with the baseline data. Key: Baseline name. Value: Dictionary with data
        id_txt:           String with dataset id.
    """
    print(f"Saving baseline files in txt/{id_txt}")

    for baseline, data in baselines.items():
        sta_1, _, sta_2 = baseline.partition("/")
        header = (
            "# "
            + f"{baseline}\n"
            + "# "
            + f"{'Date'.ljust(8)} {'Session'.ljust(7)} {'Num obs'.rjust(7)} {'Length'.rjust(14)} {'Formal error'.rjust(13)}\n"
            + "# "
            + f"{''.ljust(8)} {''.ljust(7)} {''.rjust(7)} {'[m]'.rjust(14)} {'[m]'.rjust(13)}\n"
        )
        filename = f"txt/{id_txt}/Baseline_{sta_1}_{sta_2}_{id_txt}.txt"

        dates = [d.strftime(DATE_FMT) for d in data["date"]]
        session_codes = data["session_code"]
        length = data["length"]
        ferr = data["ferr"]
        num_obs = data["num_obs"]
        
        with open(filename, "wt") as fid:
            fid.write(header)
            for d, s, n, l, e in zip(dates, session_codes, num_obs, length, ferr):
                fid.write(f"{d:10} {s:7} {n:7d} {l:14.4f} {e:13.4f}\n")


def plot_baselines(baselines, id_txt):
    """Plot baseline length for each baseline

    Plots the time evolution of the baseline length with formal errors. Each baseline is saved to
    individual plots in the 'img' directory.

    Args:
        baselines:        Dictionary with the baseline data. Key: Baseline name. Value: Dictionary with data
        id_txt:           String with dataset id.
    """
    print(f"Creating baseline plots in img/{id_txt}")
    
    for baseline, data in baselines.items():
        sta_1, _, sta_2 = baseline.partition("/")
        dates = data["date"]
        baseline_length = data["length"]
        baseline_length_ferr = data["ferr"]
        color = data["num_obs"]
        repeatability = data.get("repeatability", np.nan)
        trend = data.get("trend", [np.nan]*len(dates))
        
        fig = plt.figure(figsize=(12, 10), dpi=150)
        plt.errorbar(dates, baseline_length, yerr=baseline_length_ferr, fmt="o", marker=None, zorder=0, mew=0, ecolor="tab:gray")
        im = plt.scatter(dates, baseline_length, c=color, zorder=100)
        plt.plot(dates, trend)
        plt.gca().get_xaxis().set_major_formatter(mdates.DateFormatter(DATE_FMT))
        plt.gca().xaxis.set_major_locator(mt.LinearLocator(numticks=7))
        plt.gca().ticklabel_format(axis="y", style="plain")
        plt.gca().yaxis.get_major_formatter().set_useOffset(False)
        plt.grid(axis="y", linestyle="--")
        plt.ylabel("Baseline length [m]")
        plt.title(f"{baseline} Baseline Length Repeatability: {repeatability:6.4f} [m]")
        cbar = fig.colorbar(im, use_gridspec=True)
        cbar.set_label("Number of observations on baseline used in solution")
        plt.savefig(f"img/{id_txt}/Baseline_{sta_1}_{sta_2}_{id_txt}.png", bbox_inches="tight")
        plt.close()


def save_repeatability(baselines, id_txt="default"):
    filename = f"txt/{id_txt}/Baseline_length_repeatability_{id_txt}.txt"
    
    print(f"Saving baseline length repeatability to {filename}")
    
    with open(filename, "wt") as fid:
        header = "# " + f"{'Baseline'.ljust(18)} {'Num obs'.rjust(7)} {'Mean Length'.rjust(14)}  {'Repeatability'.rjust(13)} \n" + \
                 "# " + f"{''.ljust(18)} {''.rjust(7)} {'[m]'.rjust(14)}  {'[m]'.rjust(13)} \n"
        fid.write(header)
        for baseline, data in baselines.items():
            mean_length = data.get("mean_length")
            repeatability = data.get("repeatability")
            num_obs = data.get("repeatability_num_obs")
            if mean_length and repeatability:
                fid.write(f"{baseline:20} {num_obs:7d} {mean_length:14.4f}  {repeatability:13.4f}\n")

def plot_repeatability(baselines, id_txt="default", interactive=False):
    print(f"Creating baseline length repeatability plot in img/{id_txt}")
    
    length = []
    repeatability = []
    num_obs = []
    plotted_baselines = []
    for baseline, data in baselines.items():
        if "repeatability" in data and "mean_length" in data:
            length.append(data["mean_length"])
            repeatability.append(data["repeatability"])
            num_obs.append(data["repeatability_num_obs"])
            plotted_baselines.append(baseline)
    
    length = np.array(length) * Unit.m2km
    fig = plt.figure(figsize=(12,10), dpi=150)
    im = plt.scatter(length, repeatability, c=num_obs, picker=True)
    plt.xlabel("Baseline length [km]")
    plt.ylabel("Baseline length repeatability [m]")
    plt.grid(axis="y", linestyle="--")
    cbar = fig.colorbar(im, use_gridspec=True)
    cbar.set_label("Number of sessions with baseline used in solution")

    def print_bl(event):
        idx = event.ind[0]
        print(f"Baseline: {plotted_baselines[idx]}, Repeatbility: {repeatability[idx]} [m], Length: {length[idx]} [km], Num sessions: {num_obs[idx]}")

    fig.canvas.mpl_connect("pick_event", print_bl)
    plt.savefig(f"img/{id_txt}/Baseline_length_repeatability_{id_txt}.png", bbox_inches="tight")
    if interactive:
        plt.show()
    plt.close()
    
    # Maptest 
    plt.figure(figsize=(12,8), dpi=150)
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.set_global()
    ax.coastlines()
    cmap = plt.colormaps.get_cmap(config.there.colormap.str)
    lat = []
    lon = []
    
    
    repeat = [bl_data["repeatability"] for bl, bl_data in baselines.items() if "repeatability" in bl_data ]
    norm = colors.Normalize(vmin=min(repeat), vmax=max(repeat))
    for baseline, bl_data in baselines.items():
        if not "repeatability" in bl_data:
            continue
        sta_1, _, sta_2 = baseline.partition("/")
        sta_1_lon = np.mean(bl_data["lon_1"])
        sta_2_lon = np.mean(bl_data["lon_2"])
        sta_1_lat = np.mean(bl_data["lat_1"])
        sta_2_lat = np.mean(bl_data["lat_2"])
        repeatability = bl_data["repeatability"]
        color = cmap(norm(repeatability))
        plt.plot(
            [sta_1_lon, sta_2_lon],
            [sta_1_lat, sta_2_lat],
            c=color,
            linestyle="-",
            #linewidth=b_width[baseline],
            transform=ccrs.Geodetic(),
        )
        lon.append(sta_1_lon)
        lon.append(sta_2_lon)
        lat.append(sta_1_lat)
        lat.append(sta_2_lat)
        repeat.append(repeatability)
    plt.scatter(lon, lat, transform=ccrs.PlateCarree(), color="black", marker=".")
    plt.title(f"Baseline length repeatabilities")
    cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=plt.gca(), use_gridspec=True)
    cbar.set_label("Repeatability [m]")
    plt.savefig(f"img/{id_txt}/Baseline_length_repeatability_map{id_txt}.png", bbox_inches="tight")
    if interactive:
        plt.show()
    plt.close()
    
def compute_repeatability(baselines):
    print(f"Computing baseline length repeatability")
    
    for baseline, data in baselines.items():
        num_obs = len(data["length"])
        if num_obs > 3:
            t = Time(data["date"], scale="utc", fmt="datetime")
            l = data["length"]
            e = data["ferr"]
            #repeatability, length = wblr(l, t.mjd, e)
            repeatability, length, trend, keep_idx  = blr(l, t.mjd, e)
            data["repeatability"] = repeatability
            data["mean_length"] = length
            data["repeatability_num_obs"] = np.sum(keep_idx)
            data["trend"] = trend

            # update original data based on keep_idx
            if data["repeatability_num_obs"] < num_obs:
                print(f"Discarded baseline lengths for computation of repeatability for {baseline} for {t.date[~keep_idx]} because ferr is zero")        

        

def wblr(bl, mjd, ferr):
    """Computes weighted baseline length repeatability for short baselines

    The weighted mean of baselines are used in the the computation. Only makes
    sense for baselines lengths that do not change in the given time period. 

    Based on equation (5.10) in {hofmeister2016}

    Args:
        bl:     estimated baseline lengths for one baseline [m]
        mjd:    modified julian date for estimated baseline length
        w:      formal error of estimated baseline lengths [m]

    Returns: 
        weighted baseline length repeatability
        weighted mean baseline length
    """
    
    # TODO: discontinuities, itrf snx?

    w = 1 / ferr ** 2
    w_sum = np.sum(w)

    bl_wmean = np.sum(w * bl) / w_sum
    wblr = np.sqrt(np.sum(w * (bl - bl_wmean) ** 2) / (w_sum - np.sum(w ** 2) / w_sum))
    return wblr, bl_wmean

def blr(bl, mjd, ferr):
    """Computes baseline length repeatability for long baselines

    The mean baseline length after removing a linear trend is used in the computation. 
    Does not account for jumps in the baseline length for the given time period.
    
    Based on equation (6.3) in {hofmeister2016}

    Args:
        bl:     estimated baseline lengths for one baseline [m]
        mjd:    modified julian date for estimated baseline length
        w:      formal error of estimated baseline lengths [m]

    Returns: 
        baseline length repeatability
        mean baseline length (after detrending)
        trend
        
    """
    keep_idx = ferr != 0

    w = 1 / ferr[keep_idx] ** 2
    w_sum = np.sum(w)

    lin_func = np.polynomial.Polynomial.fit(mjd[keep_idx], bl[keep_idx], 1, w=w)
    trend = lin_func(mjd)

    bl_mean = np.mean(trend)
    blr = np.sqrt(np.sum(w * (bl[keep_idx] - trend[keep_idx]) ** 2) / (w_sum - np.sum(w ** 2) / w_sum))
    
    return blr, bl_mean, trend, keep_idx

def parse_args():
    """Define and parse input arguments to the script"""

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--id", help="Dataset id", type=str, default="")
    parser.add_argument(
        "--save_baselines",
        action="store_true",
        help="Enable this option to save all baselines results in one file per baseline",
    )
    parser.add_argument(
        "--plot_baselines",
        action="store_true",
        help="Enable this option to plot baseline vs time with one baseline per plot",
    )
    parser.add_argument(
        "--read_txt",
        action="store_true",
        help="Enable this option to read baseline results from 'txt' directory instead of session work directory",
    )
    parser.add_argument(
        "--save_repeatability",
        action="store_true",
        help="Enable this option to save the computed repeatabilities in a text file",
    )
    parser.add_argument(
        "--plot_repeatability",
        action="store_true",
        help="Enable this option to plot repeatability vs length",
    )
    parser.add_argument(
        "--delete_plots",
        action="store_true",
        help="Enable this option to delete all plots for a specific id",
        )
    parser.add_argument(
        "--delete_txt",
        action="store_true",
        help="Enable this option to delete all txt files for a specific id",
        )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Enable this option to display a interactive plot (for the plots that have this ability)"
        )
    return parser.parse_args()

def delete_files(extension, id_txt):
    print(f"Deleting *{extension} files")
    files = glob.glob(f"*/{id_txt}/*{extension}")
    for f in files:
        os.remove(f)

def main():
    # Setup config
    pipeline = "vlbi"
    config.read_pipeline(pipeline)
    config.files.profiles = [pipeline]

    args = parse_args()

    id_txt = args.id if args.id else "default"
    os.makedirs(f"img/{id_txt}/", exist_ok=True)
    os.makedirs(f"txt/{id_txt}/", exist_ok=True)

    if args.delete_plots:
        delete_files(".png", id_txt)

    if args.delete_txt:
        delete_files(".txt", id_txt)

    if args.read_txt:
        baselines = read_baselines(id_txt)
    else:
        baselines = read_session_baseline(dset_id=args.id)

    compute_repeatability(baselines)
    
    if args.plot_repeatability:
        plot_repeatability(baselines, id_txt, args.interactive)
    if args.save_repeatability:
        save_repeatability(baselines, id_txt)
    if args.plot_baselines:
        plot_baselines(baselines, id_txt)
    if args.save_baselines:
        save_baselines(baselines, id_txt)


    # Compute repeatability per baseline
    # Remove outliers?
    # Plot repeatability vs length
    

if __name__ == "__main__":
    main()
