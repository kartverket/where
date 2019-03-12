"""Transformations

Description:
------------

Transforms between GCRS and ITRS according to IERS 2010 Conventions.

References:
-----------

"""
# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import exceptions as mg_exceptions

# Where imports
from where.lib import rotation

def g2t_pos(gcrs: "GcrsPosition", time: "Time" = None) -> "TrsPosition":
    """Transforms input array from gcrs to trs coordinates"""
    if time is None:
        time = gcrs.time
        if time is None:
            raise mg_exceptions.InitializationError("Time is not defined")
    gcrs = np.asarray(gcrs)[:, :, None]
    return (rotation.gcrs2trs(time) @ gcrs)[:, :, 0]


def t2g_pos(trs: "TrsPosition", time: "Time" = None) -> "GcrsPosition":
    """Transforms input array from trs to gcrs coordinates"""
    if time is None:
        time = trs.time
        if time is None:
            raise mg_exceptions.InitializationError("Time is not defined")
    trs = np.asarray(trs)[:, :, None]
    return (rotation.trs2gcrs(time) @ trs)[:, :, 0]


def g2t_vel(gcrs: "GcrsVelocity", time: "Time" = None) -> "TrsVelocity":
    """Transforms input array from gcrs to trs coordinates"""
    if time is None:
        time = gcrs.time
        if time is None:
            raise mg_exceptions.InitializationError("Time is not defined")
    gcrs = np.asarray(gcrs)[:, :, None]
    return (rotation.dgcrs2trs_dt(time) @ gcrs)[:, :, 0]

def t2g_vel(trs: "TrsVelocity", time: "Time" = None) -> "GcrsVelocity":
    """Transforms input array from trs to gcrs coordinates"""
    if time is None:
        time = trs.time
        if time is None:
            raise mg_exceptions.InitializationError("Time is not defined")
    trs = np.asarray(trs)[:, :, None]
    return (rotation.dtrs2gcrs_dt(time) @ trs)[:, :, 0]

def g2t_posvel(gcrs: "GcrsPosVel", time: "Time" = None) -> "TrsPosVel":
    """Transforms input array from gcrs to trs coordinates"""
    if time is None:
        time = gcrs.time
        if time is None:
            raise mg_exceptions.InitializationError("Time is not defined")
    gcrs = np.asarray(gcrs)[:, :, None]
    g2t = rotation.gcrs2trs(time)
    dg2t = rotation.dgcrs2trs_dt(time)

    # Form block diagonal matrix with shape (num_obs, 6, 6)
    transformation = np.block([[g2t, np.zeros(g2t.shape)], [np.zeros(dg2t.shape), dg2t]])
    return (transformation @ gcrs)[:, :, 0]

def t2g_posvel(trs: "TrsPosVel", time: "Time" = None) -> "GcrsPosVel":
    """Transforms input array from trs to gcrs coordinates"""
    if time is None:
        time = trs.time
        if time is None:
            raise mg_exceptions.InitializationError("Time is not defined")
    trs = np.asarray(trs)[:, :, None]
    t2g = rotation.trs2gcrs(time)
    dt2g = rotation.dtrs2gcrs_dt(time)

    # Form block diagonal matrix with shape (num_obs, 6, 6)
    transformation = np.block([[t2g, np.zeros(t2g.shape)], [np.zeros(dt2g.shape), dt2g]])
    return (transformation @ trs)[:, :, 0]
