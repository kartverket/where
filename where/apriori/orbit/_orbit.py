"""An abstract baseclass for apriori orbits

Description:
------------

Implementations of apriori orbits should inherit from `AprioriOrbit` and define their own read and calculate methods.


"""
# Standard library imports
import datetime
from typing import Any, Dict

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import exceptions
from midgard.math.constant import constant

# Where imports
from where.lib import config
from where.data import dataset3 as dataset
from where.lib import log
from where.lib import util


class AprioriOrbit:
    """An abstract baseclass for initial orbits
    """

    name = "Overwritten by subclasses"

    def __init__(
        self,
        rundate: datetime.date,
        **kwargs: Dict[str, Any],
    ) -> None:
        """Set up a new AprioriOrbit object, does not parse any data

        TODO: Remove dependency on rundate, use time to read correct files. (What to do with dataset?)

        Args:
            rundate:    Date of model run.
        """
        # MURKS: Should it be done like that. The technique is normally not given for unittest routines (like
        #       test_broadcast.py).
        try:
            pipeline = config.analysis.pipeline.str
        except (AttributeError, exceptions.MissingSectionError):
            pipeline = None

        self.rundate = rundate
        self.pipeline = pipeline

        self._dset_raw = dataset.Dataset(rundate=rundate, pipeline=pipeline, stage=self.name, label="raw")
        self._dset_edit = dataset.Dataset(rundate=rundate, pipeline=pipeline, stage=self.name, label="edit")
        self._dset = None

    @property
    def dset_raw(self) -> "Dataset":
        """Dataset representing raw data from apriori orbit files

        Reads data if the raw data are not already present.
        """
        if not self._dset_raw.num_obs:
            self._read(self._dset_raw)

        return self._dset_raw

    @property
    def dset_edit(self) -> "Dataset":
        """Dataset representing edit data from apriori orbit files

        Edits data if the edit data are not already present.
        """
        if not self._dset_edit.num_obs:
            self._edit(self._dset_edit)

        return self._dset_edit

    @property
    def dset(self) -> "Dataset":
        """Dataset representing calculated apriori orbit

        Calculates data from `dset_edit` if the data are not already present.
        """
        if self._dset == None:
            self._dset = dataset.Dataset(rundate=self.rundate, pipeline=self.pipeline, stage=self.name, label="orbit")
            self._calculate(self._dset, self._dset_edit)

        return self._dset

    def calculate_orbit(self, dset: "Dataset", time: str = "time") -> None:
        """Set Dataset representing calculated apriori orbit

        Args:
            dset:   A dataset containing the data.
            time:   Define time fields to be used. It can be for example 'time' or 'sat_time'. 'time' is related to 
                    observation time and 'sat_time' to satellite transmission time.
        """
        if not dset.num_obs:
            log.fatal(f"Dataset is empty. No observation epochs given for calculating orbits.")

        self._dset = dataset.Dataset(rundate=self.rundate, pipeline=self.pipeline, stage=self.name, label="orbit")
        self._calculate(self._dset, dset, time=time)

    #
    # Abstract methods
    #
    def _read(self, dset_raw: "Dataset"):
        """Read raw data
        """
        util.not_implemented()

    def _edit(self, dset_edit: "Dataset"):
        """Edit raw data
        """
        util.not_implemented()

    def _calculate(self, dset: "Dataset"):
        """Calculate orbit data
        """
        util.not_implemented()

    #
    # Common methods for all apriori orbits
    #
    def relativistic_clock_correction(self, sat_pos: np.ndarray, sat_vel: np.ndarray) -> np.ndarray:
        """Determine relativistic clock correction due to orbit eccentricity

        The correction is caluclated for precise and broadcast orbits after Eq. 10.10 and 10.11 in :cite:`iers2010`.

        TODO: This routine should be placed in Midgard, e.g. models/ or math/?

        Args:
            sat_pos:  Array with satellite positions.
            sat_vel:  Array with satellite velocities.

        Returns:
            Relativistic clock correction due to orbit eccentricity corrections for each observation
        """
        return -2.0 / constant.c * np.einsum("ij,ij->i", sat_pos, sat_vel)
        # return -2 / constant.c * (sat_pos[:, None, :] @
        #                           sat_vel[:, :, None])[:, 0, 0]

    def satellite_clock_correction_com(
        self, antex: "AntennaCorrection", sys_freq: Dict[str, Dict[str, str]]
    ) -> np.ndarray:
        """Determine satellite clock correction related to center of mass (CoM)

        The satellite clock correction is based on Section 20.3.3.3.3.1 in :cite:`is-gps-200h`.

        Args:
            antex:     Antenna correction object based including ANTEX file data
            sys_freq:  Dictionary with frequency or frequency combination given for GNSS identifier:
                         sys_freq = { <sys_id>: <freq> }  (e.g. sys_freq = {'E': 'E1',  'G': 'L1_L2'} )

        Returns:
            GNSS satellite clock corrections for each observation in [m] related to CoM (Note: without relativistic 
            orbit eccentricity correction)
        """
        correction = antex.satellite_phase_center_offset(self.dset, sys_freq)
        return self.satellite_clock_correction(self.dset) + correction.yaw.z
