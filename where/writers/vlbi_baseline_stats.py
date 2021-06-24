"""Write statistics about baselines to file

Description:
------------

Write simple statistics about a VLBI analysis to screen and file. For each baseline, each station
and for the whole dataset the following information is printed:

* Number of observations
* Bias (average residual) in meters
* RMS of residual in meters
* Stars, ``*``, indicating high RMS'es

Example:
--------

.. highlight:: none

Example:

  Station   Station   Num obs      Offset         RMS Outlier
                            #         [m]         [m]
  ALL                    1921   -0.000006    0.008974 
  HOBART26                154    0.000257    0.014238 
  KOKEE                   631   -0.000069    0.011818 
  NYALES20                989    0.000346    0.008244 
  SVETLOE                 972   -0.000277    0.007016 
  TIGOCONC                186    0.000003    0.011267 
  WETTZELL                910   -0.000099    0.007514 
  HOBART26  KOKEE          78    0.000343    0.013261 
  HOBART26  NYALES20       19    0.001131    0.016477 
  HOBART26  SVETLOE        26   -0.003174    0.012600 
  HOBART26  TIGOCONC       22    0.000611    0.016259 
  HOBART26  WETTZELL        9    0.006708    0.016331 
  KOKEE     HOBART26       78    0.000343    0.013261 
  KOKEE     NYALES20      191    0.000900    0.012089 
  KOKEE     SVETLOE       159   -0.000502    0.011344 
  KOKEE     TIGOCONC       73   -0.000795    0.010306 
  ...

"""

# Standard library imports
import itertools
import math

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where.lib import config
from where.lib import log

MAX_NUM_STARS = 15
STAR_THRESHOLD = 0.02


@plugins.register
def baseline_stats(dset):
    """Write statistics about baselines to file.

    Args:
        dset:   Dataset, information about model run.
    """
    stats_str = ["Statistics about stations and baselines"]
    stats_str.append(f"{'Station'.ljust(9)} {'Station'.ljust(9)} {'Num obs'.rjust(7)} {'Offset'.rjust(11)} {'RMS'.rjust(11)} Outlier")
    stats_str.append(f"{''.ljust(9)} {''.ljust(9)} {'#'.rjust(7)} {'[m]'.rjust(11)} {'[m]'.rjust(11)}")
    baselines = itertools.permutations(dset.unique("station"), 2)
    idx = np.ones(dset.num_obs, dtype=bool)
    stats_str.append(_write_line("ALL", "", dset, idx))

    for sta in dset.unique("station"):
        idx = dset.filter(station=sta)
        stats_str.append(_write_line(sta, "", dset, idx))

    for sta_1, sta_2 in baselines:
        idx = np.logical_and(dset.filter(station=sta_1), dset.filter(station=sta_2))
        stats_str.append(_write_line(sta_1, sta_2, dset, idx))

    with config.files.open("output_baseline_stats", file_vars=dset.vars, mode="wt") as fid:
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

    return f"{sta_1:<9s} {sta_2:<9s} {np.sum(idx):>7d} {bias:>11.6f} {rms:>11.6f} {stars}"
