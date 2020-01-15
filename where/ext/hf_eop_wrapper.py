#!/usr/bin/env python3
"""
Where wrapper for HF EOP functions

Description:
------------

    The fortran function calc_hf_eop_xyu returns a 4x2 matrix with x, y, ut1, lod and their derivatives.
    This simple wrapper is splitting this into two functions.


References:
-----------

.. [1] Notes from the high frequency EOP working group: 
       https://ivscc.gsfc.nasa.gov/hfeop_wg/software/

.. [2] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
       IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

"""
# Standard library imports
from functools import lru_cache

# Third party imports
import numpy as np

# Midgard imports
from midgard.math.unit import Unit

# Where imports
from where.data.time import delta_ut1_utc
from where.ext import hf_eop
from where.lib import config


download_missing = config.where.files.download_missing.bool
hf_eop.hfeop_xyu.import_tides_xyu(str(config.files.path("eop_tides_desai", download_missing=download_missing)))


@lru_cache()
def hf_eop_xyu(time):
    """Computes the tidal variations of the EOP parameters x, y, ut1, lod

    Args:
            time:    epochs (see where.data.time for more info)

    Returns:
        tidal variations for polar motion x ( microarcseconds)
        tidal variations for polar motion y ( microarcseconds)
        tidal variations for ut1 (microseconds)
        tidal variations for lod (microseconds ??)
    """
    timediff = _delta_tt_ut1(time)

    if time.size == 1:
        return hf_eop.hfeop_xyu.calc_hf_eop_xyu(time.tt.mjd, timediff)[:, 0]
    else:
        # Only loop over unique epochs
        _, idx, r_idx = np.unique(np.asarray(time), return_index=True, return_inverse=True)
        return np.array(
            [hf_eop.hfeop_xyu.calc_hf_eop_xyu(t.tt.mjd, dt)[:, 0] for t, dt in zip(time[idx], timediff[idx])]
        )[r_idx]


@lru_cache()
def hf_eop_xyu_derivative(time):
    """Computes the derivative of tidal variations of the EOP parameters x, y, ut1, lod

    Args:
            time:    epochs (see where.data.time for more info)

    Returns:
        derivative of tidal variations for polar motion x ( microarcseconds / ?)
        derivative of tidal variations for polar motion y ( microarcseconds / ? )
        derivative of tidal variations for ut1 (microseconds / ? )
        derivative of tidal variations for lod (microseconds?? / ?)
    """

    timediff = _delta_tt_ut1(time)

    if time.size == 1:
        return hf_eop.hfeop_xyu.calc_hf_eop_xyu(time.tt.mjd, timediff)[:, 1]
    else:
        # Only loop over unique epochs
        _, idx, r_idx = np.unique(np.asarray(time), return_index=True, return_inverse=True)
        return np.array(
            [hf_eop.hfeop_xyu.calc_hf_eop_xyu(t.tt.mjd, dt)[:, 1] for t, dt in zip(time[idx], timediff[idx])]
        )[r_idx]


def _delta_tt_ut1(time):
    """Computes TT - UT1 in seconds"""

    # Get UT1 - UTC without applying any models to the a priori time series (to avoid infinite recursion)
    ut1_utc = delta_ut1_utc(time.utc, models=())

    # Compute TT - UTC, then substract UT1 - UTC
    jd_int_diff = (time.tt.jd_int - time.utc.jd_int) * Unit.day2second
    jd_frac_diff = (time.tt.jd_frac - time.utc.jd_frac - ut1_utc) * Unit.day2second
    return jd_int_diff + jd_frac_diff
