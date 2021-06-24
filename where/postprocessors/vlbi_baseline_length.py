"""Add estimated baseline lengths to dataset

Description:
------------

"""
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
    trf = apriori.get("trf", time=dset.time.utc.mean, reference_frames="itrf:2014, vtrf, custom, vlbi_obs")
    names = apriori.get("vlbi_station_codes")

    baselines = dset.unique("baseline")
    lengths = []
    ferrs = []
    
    save_to_file = config.tech.vlbi_baseline_length.save_to_file.bool
    bl_str = []

    
    log_and_write("# " + f"{'Baseline'.ljust(18)} {'Length'.rjust(14)} {'Formal error'.rjust(13)}", bl_str)
    log_and_write("# " + f"{''.ljust(18)} {'[m]'.rjust(14)} {'[m]'.rjust(13)}", bl_str)
    for i, bl in enumerate(baselines):
        sta_1, _, sta_2 = bl.partition("/")
        pos_apriori_1 = trf[names[sta_1]["cdp"]].pos
        pos_apriori_2 = trf[names[sta_2]["cdp"]].pos
        
        try:
            corr_1 = [
                np.mean(dset.state[f"vlbi_site_pos-{sta_1}_x"]),
                np.mean(dset.state[f"vlbi_site_pos-{sta_1}_y"]),
                np.mean(dset.state[f"vlbi_site_pos-{sta_1}_z"]),
            ]
        except AttributeError:
            log_and_write("# " + f"{bl:18} {sta_1} position not estimated", bl_str)
            continue
        
        try:
            corr_2 = [
                np.mean(dset.state[f"vlbi_site_pos-{sta_2}_x"]),
                np.mean(dset.state[f"vlbi_site_pos-{sta_2}_y"]),
                np.mean(dset.state[f"vlbi_site_pos-{sta_2}_z"]),
            ]
        except AttributeError:
            log_and_write("# " + f"{bl:18} {sta_2} position not estimated", bl_str)
            continue
        
        pos_correction_1 = PositionDelta(corr_1, system="trs", ref_pos=pos_apriori_1)
        pos_correction_2 = PositionDelta(corr_2, system="trs", ref_pos=pos_apriori_2)

        pos_1 = pos_apriori_1 + pos_correction_1
        pos_2 = pos_apriori_2 + pos_correction_2

        baseline = pos_2 - pos_1

        A = np.array(
            [
                -(pos_2.x - pos_1.x) / baseline.length,
                -(pos_2.y - pos_1.y) / baseline.length,
                -(pos_2.z - pos_1.z) / baseline.length,
                +(pos_2.x - pos_1.x) / baseline.length,
                +(pos_2.y - pos_1.y) / baseline.length,
                +(pos_2.z - pos_1.z) / baseline.length,
            ]
        )
        Q_session = dset.meta["normal equation"]["covariance"]
        session_params = dset.meta["normal equation"]["names"]
        Q_sta = extract_sta_covariance(Q_session, session_params, sta_1, sta_2)
        
        bl_length_ferr = np.sqrt(A @ Q_sta @ A.T)
        
        dset.meta.add("baseline_length", baseline.length, section=bl)
        dset.meta.add("baseline_length_ferr", bl_length_ferr, section=bl)
        dset.meta.add("__unit__", "meter", section=bl)
        
        log_and_write(f"{bl:20} {baseline.length:14.4f} {bl_length_ferr:13.4f}", bl_str)
        
        lengths.append(baseline.length)
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
    param_idx1 = np.char.startswith(param_names, f"vlbi_site_pos-{sta_1}")
    param_idx2 = np.char.startswith(param_names, f"vlbi_site_pos-{sta_2}")
    if not param_idx1.any() or not param_idx2.any():
        return None
    Q_sta = np.empty(shape=(6, 6))
    Q_sta[0:3, 0:3] = np.array(covariance)[param_idx1][:, param_idx1]
    Q_sta[0:3, 3:6] = np.array(covariance)[param_idx1][:, param_idx2]
    Q_sta[3:6, 0:3] = np.array(covariance)[param_idx2][:, param_idx1]
    Q_sta[3:6, 3:6] = np.array(covariance)[param_idx2][:, param_idx2]
    return Q_sta    
    
        

            
