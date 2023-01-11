"""Time data

Extends the Midgard time classes with additional time scales

"""
# Standard library imports
from typing import Any, Callable, Dict, List, Tuple, TypeVar

# Third party imports
import numpy as np

# Midgard imports
from midgard.data import _time as mg_time
from midgard.data.time import Time, TimeDelta, is_time, is_timedelta  # noqa
from midgard.math.constant import constant
from midgard.math.unit import Unit

# Where imports
from where import apriori
from where.ext import iers_2010 as iers
from where.ext import sofa_wrapper as sofa

# Type specification: scalar float or numpy array
np_float = TypeVar("np_float", float, np.ndarray)


######################################################################################################################
# Transitions between time scales
######################################################################################################################

def delta_ut1_utc(time: "TimeArray", models=None) -> "np_float":

    if time.scale == "ut1":
        # pretend that ut1 is utc to get an approximate delta value
        eop = apriori.get("eop", time=Time(time.mjd, fmt="mjd", scale="utc"), models=models)
        tmp_utc_mjd = time.mjd - eop.ut1_utc * Unit.second2day

        # get a better delta value with a new computed utc value
        eop = apriori.get("eop", time=Time(tmp_utc_mjd, fmt="mjd", scale="utc"), models=models)
        return -eop.ut1_utc * Unit.second2day
    else:
        # time scale is utc
        eop = apriori.get("eop", time=time, models=models)
        return eop.ut1_utc * Unit.second2day


def delta_tdb_tcb(time: "TimeArray") -> "np_float":
    if time.scale == "tdb":
        # Equation (10.3) in IERS Conventions 2010. Separate terms to avoid loss of precision
        dt = constant.T_0_jd1 - time.jd1 + constant.T_0_jd2 - time.jd2
        return constant.TDB_0 * Unit.second2day + dt * constant.L_B / (1 - constant.L_B)
    else:
        # time scale is tcb
        # Equation (10.3) in IERS Conventions 2010. Separate terms to avoid loss of precision
        dt = time.jd1 - constant.T_0_jd1 + time.jd2 - constant.T_0_jd2
        return constant.TDB_0 * Unit.second2day - dt * constant.L_B


def delta_tcb_tcg(time: "TimeArray") -> "np_float":
    # See Note 1) in hf2002_iers for time scales explanation
    if time.scale == "tcb":
        if time.size == 1:
            tcb_tcg = iers.hf2002_iers(time.tdb.jd)
        else:
            tcb_tcg = np.array([iers.hf2002_iers(t) for t in time.tdb.jd])
        return -tcb_tcg * Unit.second2day
    else:
        # time scale is tcg
        if time.size == 1:
            tcb_tcg = iers.hf2002_iers(time.tt.jd)
        else:
            tcb_tcg = np.array([iers.hf2002_iers(t) for t in time.tt.jd])
        return tcb_tcg * Unit.second2day


#
# Time scale conversions
#


def _utc2ut1(utc: "TimeArray") -> ("np_float", "np_float"):
    """Convert UTC to UT1"""
    return utc.jd1, utc.jd2 + delta_ut1_utc(utc)


def _ut12utc(ut1: "TimeArray") -> ("np_float", "np_float"):
    """Convert UT1 to UTC"""
    return ut1.jd1, ut1.jd2 + delta_ut1_utc(ut1)


def _tdb2tcb(tdb: "TimeArray") -> ("np_float", "np_float"):
    """Convert TDB to TCB"""
    return tdb.jd1, tdb.jd2 + delta_tdb_tcb(tdb)


def _tcb2tdb(tcb: "TimeArray") -> ("np_float", "np_float"):
    """Convert TCB to TDB"""
    return tcb.jd1, tcb.jd2 + delta_tdb_tcb(tcb)


def _tcg2tcb(tcg: "TimeArray") -> ("np_float", "np_float"):
    """Convert TCG to TCB"""
    return tcg.jd1, tcg.jd2 + delta_tcb_tcg(tcg)


def _tcb2tcg(tcb: "TimeArray") -> ("np_float", "np_float"):
    """Convert TCG to TT"""
    return tcb.jd1, tcb.jd2 + delta_tcb_tcg(tcb)


#
# New time scales
#
@mg_time.register_scale(convert_to=dict(utc=_ut12utc), convert_from=dict(utc=_utc2ut1))
class Ut1Time(mg_time.TimeArray):

    scale = "ut1"

    @property
    @Unit.register(("radians",), module=mg_time)
    def gmst(self):
        """Greenwich mean time based on the IAU2006 precession model"""
        return sofa.vectorized_gmst06(self)

    @property
    @Unit.register(("radians",), module=mg_time)
    def gst(self):
        """Greenwich apparent time based on the IAU2006A precession and nutation model"""
        return sofa.vectorized_gst06(self)


@mg_time.register_scale(convert_to=dict())
class Ut1TimeDelta(mg_time.TimeDeltaArray):

    scale = "ut1"


@mg_time.register_scale(convert_to=dict(tcg=_tcb2tcg), convert_from=dict(tcg=_tcg2tcb))
class TcbTime(mg_time.TimeArray):

    scale = "tcb"


@mg_time.register_scale(convert_to=dict())
class TcbTimeDelta(mg_time.TimeDeltaArray):

    scale = "tcb"


@mg_time.register_scale(convert_to=dict(tcb=_tcb2tdb), convert_from=dict(tcb=_tcb2tdb))
class TdbTime(mg_time.TimeArray):

    scale = "tdb"


@mg_time.register_scale(convert_to=dict())
class TdbTimeDelta(mg_time.TimeDeltaArray):

    scale = "tdb"

# Define shorthands for available formats, scales and conversions
Time.FORMATS = list(mg_time.TimeArray._formats().keys())
TimeDelta.FORMATS = list(mg_time.TimeDeltaArray._formats().keys())
Time.SCALES = list(mg_time.TimeArray._scales().keys())
TimeDelta.SCALES = list(mg_time.TimeDeltaArray._scales().keys())
Time.CONVERSIONS = list(mg_time.TimeArray._conversions().keys())
TimeDelta.CONVERSIONS = list(mg_time.TimeDeltaArray._conversions().keys())
