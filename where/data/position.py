""" Position classes

Extends the Midgard position classes with a additional systems

"""
# Midgard imports
from midgard.data import _position as mg_position
from midgard.data.position import Position, PositionDelta, PosVel, PosVelDelta  # noqa
from midgard.data.position import is_position, is_position_delta, is_posvel, is_posvel_delta  # noqa

# Where imports
from where.lib import transformation as trans

#
# Extra attributes
#
mg_position.register_attribute(mg_position.PositionArray, "time")
mg_position.register_attribute(mg_position.PosVelArray, "time")
mg_position.register_attribute(mg_position.VelocityArray, "time")
mg_position.register_attribute(mg_position.VelocityDeltaArray, "time")
mg_position.register_attribute(mg_position.PositionDeltaArray, "time")
mg_position.register_attribute(mg_position.PosVelDeltaArray, "time")

#
# New position systems
#
@mg_position.register_system(convert_to=dict(trs=trans.g2t_pos), convert_from=dict(trs=trans.t2g_pos))
class GcrsPosition(mg_position.PositionArray):

    system = "gcrs"
    column_names = ("x", "y", "z")
    _units = ("meter", "meter", "meter")


@mg_position.register_system(convert_to=dict(trs=trans.g2t_pos), convert_from=dict(trs=trans.t2g_pos))
class GcrsPositionDelta(mg_position.PositionDeltaArray):

    system = "gcrs"
    column_names = ("x", "y", "z")
    _units = ("meter", "meter", "meter")


@mg_position.register_system(convert_to=dict(trs=trans.g2t_posvel), convert_from=dict(trs=trans.t2g_posvel))
class GcrsPosVel(mg_position.PosVelArray):

    system = "gcrs"
    column_names = ("x", "y", "z", "vx", "vy", "vz")
    _units = ("meter", "meter", "meter", "meter/second", "meter/second", "meter/second")


@mg_position.register_system(convert_to=dict(trs=trans.g2t_posvel), convert_from=dict(trs=trans.t2g_posvel))
class GcrsPosVelDelta(mg_position.PosVelDeltaArray):

    system = "gcrs"
    column_names = ("x", "y", "z", "vx", "vy", "vz")
    _units = ("meter", "meter", "meter", "meter/second", "meter/second", "meter/second")


@mg_position.register_system(convert_to=dict(trs=trans.g2t_vel), convert_from=dict(trs=trans.t2g_vel))
class GcrsVelocity(mg_position.VelocityArray):

    system = "gcrs"
    column_names = ("vx", "vy", "vz")
    _units = ("meter/second", "meter/second", "meter/second")


@mg_position.register_system(convert_to=dict(trs=trans.g2t_vel), convert_from=dict(trs=trans.t2g_vel))
class GcrsVelocityDelta(mg_position.VelocityDeltaArray):

    system = "gcrs"
    column_names = ("vx", "vy", "vz")
    _units = ("meter/second", "meter/second", "meter/second")


@mg_position.register_system(
    convert_to=dict(trs=trans.delta_y2t_posvel), convert_from=dict(trs=trans.delta_t2y_posvel)
)
class YawPosVelDelta(mg_position.PosVelDeltaArray):

    system = "yaw"
    column_names = ("x", "y", "z", "vx", "vy", "vz")
    _units = ("meter", "meter", "meter", "meter/second", "meter/second", "meter/second")


@mg_position.register_system(convert_to=dict(trs=trans.delta_y2t), convert_from=dict(trs=trans.delta_t2y))
class YawPositionDelta(mg_position.PositionDeltaArray):

    system = "yaw"
    column_names = ("x", "y", "z")
    _units = ("meter", "meter", "meter")


@mg_position.register_system(convert_to=dict(trs=trans.delta_y2t), convert_from=dict(trs=trans.delta_t2y))
class YawVelocityDelta(mg_position.VelocityDeltaArray):

    system = "yaw"
    column_names = ("vx", "vy", "vz")
    _units = ("meter/second", "meter/second", "meter/second")
