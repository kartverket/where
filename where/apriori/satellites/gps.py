"""General information about GPS satellites

Description:



"""
# Standard library imports
import math

# External library imports
import numpy as np

# Where imports
from where import apriori
from where.lib import plugins
from where.apriori.satellites import satellite


class Gps:
    """General properties of GPS satellites

    This class should be used as a mixin together with satellite.Satellite. That is, any GPS-satellite object should
    inherit from both Gps and satellite.Satellite (in that order).
    """
    area = math.pi * 0.30 ** 2
    mass = 40.965

    def initial_posvel(self, initial_dt):
        """Calculate initial position and velocity for GPS satellites
        """
        orb = apriori.get("gnss_orbit", rundate=initial_dt)[self.name]["pos"]
        pos = orb.gcrs[0]
        vel = (orb.gcrs[1] - orb.gcrs[0]) / (orb.utc_dt[1] - orb.utc_dt[0]).total_seconds()
        return np.hstack((pos, vel))


@plugins.register
class G01(Gps, satellite.Satellite):
    pass


@plugins.register
class G02(Gps, satellite.Satellite):
    pass


@plugins.register
class G03(Gps, satellite.Satellite):
    pass


@plugins.register
class G04(Gps, satellite.Satellite):
    pass


@plugins.register
class G05(Gps, satellite.Satellite):
    pass


@plugins.register
class G06(Gps, satellite.Satellite):
    pass


@plugins.register
class G07(Gps, satellite.Satellite):
    pass


@plugins.register
class G08(Gps, satellite.Satellite):
    pass


@plugins.register
class G09(Gps, satellite.Satellite):
    pass


@plugins.register
class G10(Gps, satellite.Satellite):
    pass


@plugins.register
class G11(Gps, satellite.Satellite):
    pass


@plugins.register
class G12(Gps, satellite.Satellite):
    pass


@plugins.register
class G13(Gps, satellite.Satellite):
    pass


@plugins.register
class G14(Gps, satellite.Satellite):
    pass


@plugins.register
class G15(Gps, satellite.Satellite):
    pass


@plugins.register
class G16(Gps, satellite.Satellite):
    pass


@plugins.register
class G17(Gps, satellite.Satellite):
    pass


@plugins.register
class G18(Gps, satellite.Satellite):
    pass


@plugins.register
class G19(Gps, satellite.Satellite):
    pass


@plugins.register
class G20(Gps, satellite.Satellite):
    pass


@plugins.register
class G21(Gps, satellite.Satellite):
    pass


@plugins.register
class G22(Gps, satellite.Satellite):
    pass


@plugins.register
class G23(Gps, satellite.Satellite):
    pass


@plugins.register
class G24(Gps, satellite.Satellite):
    pass


@plugins.register
class G25(Gps, satellite.Satellite):
    pass


@plugins.register
class G26(Gps, satellite.Satellite):
    pass


@plugins.register
class G27(Gps, satellite.Satellite):
    pass


@plugins.register
class G28(Gps, satellite.Satellite):
    pass


@plugins.register
class G29(Gps, satellite.Satellite):
    pass


@plugins.register
class G30(Gps, satellite.Satellite):
    pass


@plugins.register
class G31(Gps, satellite.Satellite):
    pass


@plugins.register
class G32(Gps, satellite.Satellite):
    pass
