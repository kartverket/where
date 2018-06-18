"""An abstract baseclass for apriori orbits

Description:
------------

Implementations of apriori orbits should inherit from `AprioriOrbit` and define their own read and calculate methods.




"""

# Standard library imports

# External library imports
import numpy as np

# Where imports
from where.lib import config
from where.lib import constant
from where import data
from where.lib import util


class AprioriOrbit():
    """An abstract baseclass for initial orbits
    """

    name = "Overwritten by subclasses"

    def __init__(self, rundate, time, satellite, system=None, **kwargs):
        """Set up a new AprioriOrbit object, does not parse any data

        TODO: Remove dependency on rundate, use time to read correct files. (What to do with dataset?)

        Args:
            rundate (datetime.date):  Date of model run.
            time (TimeTable):         Time epochs at the satellite for which to calculate the apriori orbit.
            satellite (list):         Strings with names of satellites.
            system (list):            Strings with codes of systems (G, E, R, etc.).
        """
        # MURKS: Should it be done like that. The technique is normally not given for unittest routines (like
        #       test_broadcast.py).
        try:
            tech = config.analysis.tech.str
        except AttributeError:
            tech = None

        self.time = time
        self.satellite = satellite
        self.system = [s[0] for s in satellite] if system is None else system
        self._dset_raw = data.Dataset(
            rundate=rundate, tech=tech, stage=self.name, dataset_name="raw", dataset_id=0, empty=True, session=""
        )
        self._dset_edit = data.Dataset(
            rundate=rundate, tech=tech, stage=self.name, dataset_name="edit", dataset_id=0, empty=True, session=""
        )
        self._dset = data.Dataset(
            rundate=rundate, tech=tech, stage=self.name, dataset_name="orbit", dataset_id=0, empty=True, session=""
        )

    @property
    def dset_raw(self):
        """Dataset representing raw data from apriori orbit files

        Reads data if the raw data are not already present.
        """
        if not self._dset_raw.num_obs:
            self._dset_raw = self._read(self._dset_raw)

        return self._dset_raw

    @property
    def dset_edit(self):
        """Dataset representing raw data from apriori orbit files

        Reads data if the raw data are not already present.
        """
        if not self._dset_edit.num_obs:
            self._dset_edit = self._edit(self._dset_edit)

        return self._dset_edit

    @property
    def dset(self):
        """Dataset representing calculated apriori orbit

        Calculates data from `dset_raw` if the data are not already present.
        """
        if not self._dset.num_obs:
            self._dset = self._calculate(self._dset)

        return self._dset

    def update_time(self, time, satellite=None):
        """Update which time epochs to calculate orbit for

        Will trigger a recalculation of the calculated apriori orbit.

        Args:
            time (TimeTable):  New time epochs at the satellite for which to calculate the apriori orbit.
            satellite (list):  Optional updated list of satellites.
        """
        self.time = time
        self._dset = data.Dataset(
            rundate=self._dset.rundate,
            tech=config.analysis.tech.str,
            stage=self.name,
            dataset_name="orbit",
            dataset_id=self._dset.dataset_id + 1,
            empty=True,
        )

        if satellite is not None:
            self.satellite = satellite

    #
    # Abstract methods
    #
    def _read(self, dset_raw):
        """Read raw data
        """
        util.not_implemented()

    def _edit(self, dset_edit):
        """Edit raw data
        """
        util.not_implemented()

    def _calculate(self, dset):
        """Calculate orbit data
        """
        util.not_implemented()

    #
    # Common methods for all apriori orbits
    #
    def relativistic_clock_correction(self):
        """Determine relativistic clock correction due to orbit eccentricity

        The correction is caluclated for precise and broadcast orbits after Eq. 10.10 and 10.11 in :cite:`iers2010`.

        Returns:
            numpy.ndarray:    Relativistic clock correction due to orbit eccentricity corrections for each observation
        """
        return -2.0 / constant.c * np.einsum("ij,ij->i", self.dset.sat_posvel.itrs_pos, self.dset.sat_posvel.itrs_vel)
        # return -2 / constant.c * (self.dset.sat_posvel.itrs_pos[:, None, :] @
        #                           self.dset.sat_posvel.itrs_vel[:, :, None])[:, 0, 0]

    def satellite_clock_correction_com(self, antex):
        """Determine satellite clock correction related to center of mass (CoM)

        The satellite clock correction is based on Section 20.3.3.3.3.1 in :cite:`is-gps-200h`.

        Returns:
            numpy.ndarray:    GNSS satellite clock corrections for each observation in [m] related to CoM 
                              (Note: without relativistic orbit eccentricity correction)
        """
        z_yaw = self.dset.sat_posvel.convert_itrs_to_yaw(antex.satellite_phase_center_offset(self.dset))[:, 2]
        return self.satellite_clock_correction() + z_yaw
