#!/usr/bin/env python3
"""
Where wrapper for IERS 2010 functions

Description:
------------

Contains vectorized wrappers for some functions in where.ext.iers_2010.

References:
-----------

.. [2] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
       IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

"""
# Standard library imports
from functools import lru_cache

# Third party imports
import numpy as np

# Where imports
from where.ext import iers_2010 as iers


@lru_cache()
def ortho_eop(time):
    """Computes the tidal variations of the EOP parameters x, y and ut1

    Args:
            time:    epochs (see where.data.time for more info)

    Returns:
        tidal variations for polar motion x ( microarcseconds)
        tidal variations for polar motion y ( microarcseconds)
        tidal variations for ut1 (microseconds)
    """
    if time.size == 1:
        return iers.ortho_eop(time.tt.mjd)
    else:
        # Only loop over unique epochs
        _, idx, r_idx = np.unique(np.asarray(time), return_index=True, return_inverse=True)
        return np.array([iers.ortho_eop(t.tt.mjd) for t in time[idx]])[r_idx]


@lru_cache()
def utlibr(time):
    """Computes sub-diurnal libration for the EOP parameters ut1 and lod

    Args:
            time:    epochs (see where.data.time for more info)

    Returns:
        libration for ut1 (microseconds / ? )
        libration for lod (microseconds per day)
    """

    if time.size == 1:
        return iers.utlibr(time.tt.mjd)
    else:
        # Only loop over unique epochs
        _, idx, r_idx = np.unique(np.asarray(time), return_index=True, return_inverse=True)
        return np.array([iers.utlibr(t.tt.mjd) for t in time[idx]])[r_idx]


@lru_cache()
def pmsdnut2(time):
    """Computes diunal libration of the EOP parameters x and y

    Args:
            time:    epochs (see where.data.time for more info)

    Returns:
        libration for polar motion x ( microarcseconds)
        libration for polar motion y ( microarcseconds)
    """
    if time.size == 1:
        return iers.pmsdnut2(time.tt.mjd)
    else:
        # Only loop over unique epochs
        _, idx, r_idx = np.unique(np.asarray(time), return_index=True, return_inverse=True)
        return np.array([iers.pmsdnut2(t.tt.mjd) for t in time[idx]])[r_idx]


@lru_cache()
def rg_zont2(time):
    """Computes tidal variations in the Earth’s rotation with periods from 5 days to 18.6 years

    Args:
            time:    epochs (see where.data.time for more info)

    Returns:
        tidal variations for ut1 (seconds)
        tidal variations for lod (seconds per day)
        tidal variations for earth rotation speed (radians per second)
    """
    # # Julian centuries since J2000
    t_julian_centuries = (time.tt.jd - 2_451_545.0) / 36525

    if time.size == 1:
        return iers.rg_zont2(t_julian_centuries)
    else:
        # Only loop over unique epochs
        _, idx, r_idx = np.unique(np.asarray(time), return_index=True, return_inverse=True)
        return np.array([iers.rg_zont2(t) for t in t_julian_centuries[idx]])[r_idx]
