"""Transformations

Description:
------------

Transforms between GCRS and ITRS according to IERS 2010 Conventions.

References:
-----------

"""
# Standard library imports
from functools import lru_cache

# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import exceptions as mg_exceptions
from midgard.math import nputil

# Where imports
from where.lib import rotation


@nputil.hashable
@lru_cache()
def _matmul(a, b):
    return np.squeeze(a @ b)


def g2t_pos(gcrs: "GcrsPosition", time: "Time" = None) -> "TrsPosition":
    """Transforms input array from gcrs to trs coordinates"""
    if time is None:
        time = gcrs.time
        if time is None:
            raise mg_exceptions.InitializationError("Time is not defined")
    gcrs = nputil.col(gcrs)
    return _matmul(rotation.gcrs2trs(time), gcrs)


def t2g_pos(trs: "TrsPosition", time: "Time" = None) -> "GcrsPosition":
    """Transforms input array from trs to gcrs coordinates"""
    if time is None:
        time = trs.time
        if time is None:
            raise mg_exceptions.InitializationError("Time is not defined")
    trs = nputil.col(trs)
    return _matmul(rotation.trs2gcrs(time), trs)


def g2t_vel(gcrs: "GcrsVelocity", time: "Time" = None) -> "TrsVelocity":
    """Transforms input array from gcrs to trs coordinates"""
    if time is None:
        time = gcrs.time
        if time is None:
            raise mg_exceptions.InitializationError("Time is not defined")
    gcrs = nputil.col(gcrs)
    return _matmul(rotation.dgcrs2trs_dt(time), gcrs)


def t2g_vel(trs: "TrsVelocity", time: "Time" = None) -> "GcrsVelocity":
    """Transforms input array from trs to gcrs coordinates"""
    if time is None:
        time = trs.time
        if time is None:
            raise mg_exceptions.InitializationError("Time is not defined")
    trs = nputil.col(trs)
    return _matmul(rotation.dtrs2gcrs_dt(time), trs)


def g2t_posvel(gcrs: "GcrsPosVel", time: "Time" = None) -> "TrsPosVel":
    """Transforms input array from gcrs to trs coordinates"""
    if time is None:
        time = gcrs.time
        if time is None:
            raise mg_exceptions.InitializationError("Time is not defined")
    gcrs = nputil.col(gcrs)
    g2t = rotation.gcrs2trs(time)
    dg2t = rotation.dgcrs2trs_dt(time)

    transformation = np.block([[g2t, np.zeros(g2t.shape)], [dg2t, g2t]])
    return _matmul(transformation, gcrs)


def t2g_posvel(trs: "TrsPosVel", time: "Time" = None) -> "GcrsPosVel":
    """Transforms input array from trs to gcrs coordinates"""
    if time is None:
        time = trs.time
        if time is None:
            raise mg_exceptions.InitializationError("Time is not defined")
    trs = nputil.col(trs)
    t2g = rotation.trs2gcrs(time)
    dt2g = rotation.dtrs2gcrs_dt(time)

    transformation = np.block([[t2g, np.zeros(t2g.shape)], [dt2g, t2g]])
    return _matmul(transformation, trs)


def delta_t2y(trs: "TrsPositionDelta") -> "YawPositionDelta":
    """Convert position deltas from TRS to YAW"""
    t2y = rotation.trs2yaw(trs.ref_pos, trs.time)
    return _matmul(t2y, trs.mat)


def delta_y2t(yaw: "YawPositionDelta") -> "TrsPositionDelta":
    """Convert position deltas from YAW to TRS"""
    y2t = rotation.yaw2trs(yaw.ref_pos, yaw.time)
    return _matmul(y2t, yaw.mat)


def delta_t2y_posvel(trs: "TrsPosVelDelta") -> "YawPosVelDelta":
    """Convert position deltas from TRS to YAW"""
    t2y = rotation.trs2yaw(trs.ref_pos, trs.time)
    # TODO: verify this tranformation
    trs2yaw = np.block([[t2y, np.zeros(t2y.shape)], [np.zeros(t2y.shape), t2y]])
    return _matmul(trs2yaw, trs.mat)


def delta_y2t_posvel(yaw: "YawPosVelDelta") -> "TrsPosVelDelta":
    """Convert position deltas from YAW to TRS"""
    y2t = rotation.yaw2trs(yaw.ref_pos, yaw.time)
    # TODO: verify this tranformation
    yaw2trs = np.block([[y2t, np.zeros(y2t.shape)], [np.zeros(y2t.shape), y2t]])
    return _matmul(yaw2trs, yaw.mat)
