"""Basic information about satellites, extended by individual satellites

Description:
------------

This is the base class defining basic properties of satellites. Any particular satellite should be defined in
its own class, inheriting from this one.



"""

# External library imports
import numpy as np

# Where imports
from where import apriori


class Satellite:
    """General properties for all satellites

    This class should be inherited by all satellites, possibly with more specialized classes (e.g. Slr, Gps) as mixins.
    """

    short_name: str = ""
    cospar_id: int = 0
    area: float = 0
    mass: float = 0

    def __init__(self):
        """Initialize satellite object

        The satellite name is by default taken to be the lower case of the class name.
        """
        self.name = self.__class__.__name__.lower()


class Slr:
    """General properties of SLR satellites

    This class should be used as a mixin together with satellite.Satellite. That is, any SLR-satellite object should
    inherit from both Glr and satellite.Satellite (in that order).
    """

    def initial_posvel(self, initial_dt):
        """Estimate initial position and velocity of SLR satellite

        Args:
            initial_dt:  Datetime for which position and velocity is found.

        Returns:
            Position and velocity of satellite as 6-vector.
        """
        predictions = apriori.get("slr_ephemeris", rundate=initial_dt, sat_name=self.name)
        fields = ["initial_pos", "initial_vel"]

        return np.hstack([predictions[f] for f in fields])
