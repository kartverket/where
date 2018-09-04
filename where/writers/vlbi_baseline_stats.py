"""Write statistics about baselines to file

Description:
------------

Write simple statistics about a VLBI analysis to screen and file. For each baseline, each station
and for the whole dataset the following information is printed:

* Number of observations
* Length of baseline in meters (only for baselines)
* Bias (average residual) in meters
* RMS of residual in meters
* Stars, ``*``, indicating high RMS'es

Example:
--------

.. highlight:: none

The following example shows output from a VLBI analysis of three stations::

    ALL                  2686         0.00   -0.00143    0.09186 *
    HOBART12             1798         0.00    0.00072    0.08387 *
    KATH12M              1782         0.00   -0.00247    0.09918 *
    YARRA12M             1792         0.00   -0.00255    0.09196 *
    HOBART12  KATH12M     894   3431879.01    0.00082    0.09166 *
    HOBART12  YARRA12M    904   3211335.64    0.00063    0.07539 *
    KATH12M   HOBART12    894   3431879.01    0.00082    0.09166 *
    KATH12M   YARRA12M    888   2360367.26   -0.00579    0.10621 **
    YARRA12M  HOBART12    904   3211335.64    0.00063    0.07539 *
    YARRA12M  KATH12M     888   2360367.26   -0.00579    0.10621 **
"""

# Standard library imports
import itertools
import math

# External library imports
import numpy as np

# Where imports
from where.lib import files
from where.lib import log
from where.lib import plugins

MAX_NUM_STARS = 15
STAR_THRESHOLD = 0.02


@plugins.register
def baseline_stats(dset):
    """Write statistics about baselines to file.

    Args:
        dset:   Dataset, information about model run.
    """
    stats_str = ["Statistics about stations and baselines"]
    baselines = itertools.permutations(dset.unique("station"), 2)
    idx = np.ones(dset.num_obs, dtype=bool)
    stats_str.append(_write_line("ALL", "", dset, idx))

    for sta in dset.unique("station"):
        idx = dset.filter(station=sta)
        stats_str.append(_write_line(sta, "", dset, idx))

    for sta_1, sta_2 in baselines:
        idx = np.logical_and(dset.filter(station=sta_1), dset.filter(station=sta_2))
        stats_str.append(_write_line(sta_1, sta_2, dset, idx))

    with files.open("output_baseline_stats", file_vars=dset.vars, mode="wt") as fid:
        fid.write("\n".join(stats_str))
    log.out("\n  ".join(stats_str))


def _write_line(sta_1, sta_2, dset, idx):
    """Write one line of information about a station or baseline.

    At the moment, we also print the line to the screen for convenience, this might need to be
    changed when we run operationally?

    Args:
        sta_1:  String, name of Station 1.
        sta_2:  String, name of Station 2.
        dset:   Dataset, data from the model run.
        idx:    Bool-array, True for observations that should be included.
        fid:    File-pointer, file to write to.

    """
    rms = dset.rms("residual", idx=idx) if any(idx) else 0
    bias = dset.mean("residual", idx=idx) if any(idx) else 0
    stars = "" if rms < STAR_THRESHOLD else "*" * min(MAX_NUM_STARS, math.ceil(math.log2(rms / STAR_THRESHOLD)))
    bl_len = np.linalg.norm(dset.site_pos_2.itrs[idx][0] - dset.site_pos_1.itrs[idx][0]) if (sta_2 and any(idx)) else 0

    return f"{sta_1:<9s} {sta_2:<9s} {np.sum(idx):>5d} {bl_len:>11.1f} {bias:>11.6f} {rms:>11.6f} {stars}"
