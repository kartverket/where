"""Add estimated baseline lengths to dataset

Description:
------------

"""
# Standard library imports
import itertools

# External library imports
import matplotlib.pyplot as plt
import numpy as np
import os

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit
from midgard.math.constant import constant

# Where imports
from where import apriori
from where.data.position import PositionDelta
from where.lib import config
from where.lib import log

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def vlbi_baseline_length(dset: "Dataset") -> None:
    """Add estimated baseline lengths to dataset

    Station positions needs to be estimated to compute baseline lengths. 

    Args:
        dset:     A Dataset containing model data.
    """
    trf = apriori.get("trf", time=dset.time.utc.mean)

    baselines = itertools.combinations(dset.unique("station"), 2)

    lengths = []
    ferrs = []
    
    save_to_file = config.tech.vlbi_baseline_length.save_to_file.bool
    bl_str = []

    log_and_write("# " + f"{'Baseline'.ljust(18)} {'Num obs'.ljust(7)} {'Length'.rjust(14)} {'Formal error'.rjust(13)} " +
                  f"{'Latitude sta_1'.rjust(15)} {'Latitude sta_2'.rjust(15)} {'Longitude sta_1'.rjust(15)} {'Longitude sta_2'.rjust(15)}",
                  bl_str)
    log_and_write("# " + f"{'sta_1/sta_2'.ljust(18)} {''.rjust(7)} {'[m]'.rjust(14)} {'[m]'.rjust(13)} " +
                  f"{'[Degrees]'.rjust(15)} {'[Degrees]'.rjust(15)} {'[Degrees]'.rjust(15)} {'[Degrees]'.rjust(15)}",
                  bl_str)

    param_groups = set([f.split("-")[0] for f in dset.state.fields])
    Q_session = dset.meta["normal equation"]["covariance"]
    neq_params = dset.meta["normal equation"]["names"]

    for i, (sta_1, sta_2) in enumerate(baselines):
        # Name the baseline in alphabetical station order to be consistent across sessions
        bl = sta_1 + "/" + sta_2
        bl_sorted = "/".join(sorted([sta_1, sta_2]))
        bl_rsorted = "/".join(sorted([sta_1, sta_2], reverse=True))
        sta_1_sorted, _, sta_2_sorted = bl_sorted.partition("/")

        pos_apriori_1 = trf[dset.meta["station"][sta_1]["cdp"]].pos
        pos_apriori_2 = trf[dset.meta["station"][sta_2]["cdp"]].pos
        lat_1 = np.degrees(dset.meta["station"][sta_1_sorted]["latitude"])
        lat_2 = np.degrees(dset.meta["station"][sta_2_sorted]["latitude"])
        lon_1 = np.degrees(dset.meta["station"][sta_1_sorted]["longitude"])
        lon_2 = np.degrees(dset.meta["station"][sta_2_sorted]["longitude"])
        
        if "vlbi_site_pos" in param_groups:
            # Compute baseline length based on estimated station positions
            try:
                x_idx = neq_params.index(f"vlbi_site_pos-{sta_1}_x")
                y_idx = neq_params.index(f"vlbi_site_pos-{sta_1}_y")
                z_idx = neq_params.index(f"vlbi_site_pos-{sta_1}_z")

                corr_1 = [
                    dset.meta["normal equation"]["solution"][x_idx],
                    dset.meta["normal equation"]["solution"][y_idx],
                    dset.meta["normal equation"]["solution"][z_idx],
                ]
            except AttributeError:
                log_and_write("# " + f"{bl_sorted:18} {sta_1} position not estimated", bl_str)
                continue

            try:
                x_idx = neq_params.index(f"vlbi_site_pos-{sta_2}_x")
                y_idx = neq_params.index(f"vlbi_site_pos-{sta_2}_y")
                z_idx = neq_params.index(f"vlbi_site_pos-{sta_2}_z")

                corr_2 = [
                    dset.meta["normal equation"]["solution"][x_idx],
                    dset.meta["normal equation"]["solution"][y_idx],
                    dset.meta["normal equation"]["solution"][z_idx],
                ]
            except AttributeError:
                log_and_write("# " + f"{bl_sorted:18} {sta_2} position not estimated", bl_str)
                continue

            pos_correction_1 = PositionDelta(corr_1, system="trs", ref_pos=pos_apriori_1)
            pos_correction_2 = PositionDelta(corr_2, system="trs", ref_pos=pos_apriori_2)

            pos_1 = pos_apriori_1 + pos_correction_1
            pos_2 = pos_apriori_2 + pos_correction_2

            baseline = pos_2 - pos_1
            bl_length = baseline.length

            A = np.array(
                [
                    -(pos_2.x - pos_1.x) / bl_length,
                    -(pos_2.y - pos_1.y) / bl_length,
                    -(pos_2.z - pos_1.z) / bl_length,
                    +(pos_2.x - pos_1.x) / bl_length,
                    +(pos_2.y - pos_1.y) / bl_length,
                    +(pos_2.z - pos_1.z) / bl_length,
                ]
            )

            Q_sta = extract_sta_covariance(Q_session, neq_params, sta_1, sta_2)
            bl_length_ferr = np.sqrt(A @ Q_sta @ A.T)

        elif "vlbi_baseline" in param_groups:
            # Compute baseline length based on estimated baselines
            neq_params = dset.meta["normal equation"]["names"]
            if f"vlbi_baseline-{bl_sorted}" in neq_params:
                param_idx = neq_params.index(f"vlbi_baseline-{bl_sorted}")
            elif f"vlbi_baseline-{bl_rsorted}" in neq_params:
                param_idx = neq_params.index(f"vlbi_baseline-{bl_rsorted}")
            else:
                log_and_write(f"# Baseline {bl_sorted} not estimated")
                continue

            apriori_bl_length = (pos_apriori_2 - pos_apriori_1).length
            # todo: display units
            bl_length = apriori_bl_length + dset.meta["normal equation"]["solution"][param_idx]
            bl_length_ferr = np.sqrt(dset.meta["normal equation"]["covariance"][param_idx][param_idx])

        else:
            log_and_write("No station position or baseline lengths are estimated")
            return

        dset.meta.add("baseline_length", bl_length, section=bl_sorted)
        dset.meta.add("baseline_length_ferr", bl_length_ferr, section=bl_sorted)
        dset.meta.add("__unit__", "meter", section=bl_sorted)
        
        num_obs = dset.num(baseline=bl_rsorted) or dset.num(baseline=bl_sorted)
        
        log_and_write(f"{bl_sorted:20} {num_obs:7d} {bl_length:14.4f} {bl_length_ferr:13.4f} " +
                      f"{lat_1:15.10f} {lat_2:15.10f} {lon_1:15.10f} {lon_2:15.10f}", bl_str)
        
        lengths.append(bl_length)
        ferrs.append(bl_length_ferr)
        
    if config.tech.vlbi_baseline_length.plot_length_vs_ferr.bool:
        plot_length_vs_ferr(dset, np.array(lengths), np.array(ferrs))

    if save_to_file:
        with config.files.open("vlbi_baseline_lengths", file_vars=dset.vars, mode="wt") as fid:
            for line in bl_str:
                fid.write(line + '\n')


def log_and_write(text, bl_str):
    log.info(text)
    bl_str.append(text)

def plot_length_vs_ferr(dset, lengths, ferrs):
    fig = plt.figure(figsize=(12,10), dpi=150)
    plt.scatter(lengths * Unit.meter2kilometer, ferrs)
    plt.xlabel(f"Estimated baseline length [km]")
    plt.ylabel(f"Formal error [m]")
    filename = config.files.path("vlbi_baseline_length_vs_ferr", file_vars=dset.vars)
    if not filename.parent.exists(): 
        os.makedirs(filename.parent)
    plt.savefig(filename, bbox_inches='tight')
    plt.close()


def extract_sta_covariance(covariance, param_names, sta_1, sta_2):
    param_idx1 = np.char.startswith(param_names, f"vlbi_site_pos-{sta_1}_")
    param_idx2 = np.char.startswith(param_names, f"vlbi_site_pos-{sta_2}_")
    if not param_idx1.any() or not param_idx2.any():
        return None
    Q_sta = np.empty(shape=(6, 6))
    Q_sta[0:3, 0:3] = np.array(covariance)[param_idx1][:, param_idx1]
    Q_sta[0:3, 3:6] = np.array(covariance)[param_idx1][:, param_idx2]
    Q_sta[3:6, 0:3] = np.array(covariance)[param_idx2][:, param_idx1]
    Q_sta[3:6, 3:6] = np.array(covariance)[param_idx2][:, param_idx2]
    return Q_sta    
    
        

            
