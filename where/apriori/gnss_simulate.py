"""Simulating GNSS data

Description:
------------

Example:
    from datetime import datetime
    from where import apriori    

    # Get simulated data object
    sim = apriori.get('gnss_simulate', rundate=dateteime(2019,11,02), sampling_rate=300, systems=['G'], satellites=['G01'])

    # Write calculated Dataset to file
    sim.dset.write()


"""
# Standard library imports
from datetime import datetime, timedelta
from typing import List, Union

# External library imports
import numpy as np

# Midgard imports
from midgard.math.constant import constant
from midgard.dev import plugins

# Where imports
from where import apriori
from where.data import dataset3 as dataset
from where.lib import config
from where.lib import util

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
class GnssSimulate:
    """Simulating GNSS data
    """

    def __init__(
        self,
        rundate: datetime,
        sampling_rate: Union[float, None] = None,
        systems: Union[List[str], None] = None,
        satellites: Union[List[str], None] = None,
        site_pos: Union[List[float], None] = None,
        station: Union[str, None] = None,
    ):
        """Set up a new GnssSimulate object

        TODO: Remove dependency on rundate, use time to read correct files. (What to do with dataset?)

        Args:
            rundate (datetime.date):  Date of model run.
            time (TimeTable):         Time epochs at the satellite for which to calculate the apriori orbit.
            satellite (list):         Strings with names of satellites.
        """
        # MURKS: Should it be done like that. The technique is normally not given for unittest routines (like
        #       test_broadcast.py).
        try:
            pipeline = config.tech.pipeline.str
        except:
            pipeline = "gnss"

        # +MURKS
        systems = ["G"]
        satellites = ["G01", "G02"]
        site_pos = [3348186.1150, 465040.8615, 5390738.0919]
        station = "krss"
        # -MURKS
        self.sampling_rate = sampling_rate if sampling_rate else config.tech.simulate.sampling_rate.float
        self.systems = systems if systems else config.tech.simulate.systems.list
        self.satellites = satellites if satellites else config.tech.simulate.satellites.list
        self.site_pos = np.array(site_pos) if site_pos else np.array(config.tech.simulate.site_pos.list)
        self.station = station if station else config.tech.simulate.station.str
        # TODO: Generate station position

        self._dset_raw = dataset.Dataset(rundate=rundate, pipeline=pipeline, stage="simulate", label="raw")
        self._dset_orbit = dataset.Dataset(rundate=rundate, pipeline=pipeline, stage="simulate", label="orbit")
        self._dset = dataset.Dataset(rundate=rundate, pipeline=pipeline, stage="simulate", label="simulate")

    @property
    def dset_raw(self):
        """Dataset representing simulated raw data
        """
        if not self._dset_raw.num_obs:
            self._setup(self._dset_raw)

        return self._dset_raw

    @property
    def dset_orbit(self):
        """Dataset representing orbit data based on simulated raw data
        """
        if not self._dset_orbit.num_obs:
            self._orbit(self._dset_orbit)

        return self._dset_orbit

    @property
    def dset(self):
        """Dataset representing calculated apriori orbit

        Calculates data from `dset_edit` if the data are not already present.
        """
        if not self._dset.num_obs:
            self._simulate(self._dset)

        return self._dset

    @property
    def dset(self):
        """Dataset representing calculated apriori orbit

        Calculates data from `dset_edit` if the data are not already present.
        """
        if not self._dset.num_obs:
            self._simulate(self._dset)

        return self._dset

    def _orbit(self, dset: "Dataset") -> None:
        """Generate orbit data based on simulated raw data

        Args:
            dset:   A dataset containing the data.
        """
        orbit = apriori.get(
            "orbit",
            rundate=self.dset_raw.analysis["rundate"],
            system=tuple(self.dset_raw.unique("system")),
            station=self.station,
        )
        # orbit.dset_raw.write_as(stage=stage, session=station, dataset_name="raw")
        # orbit.dset_edit.write_as(stage=stage, session=station, dataset_name="edit")
        orbit.calculate_orbit(self.dset_raw, time="time")
        dset = self.dset_raw.copy()

    def _setup(self, dset: "Dataset") -> None:
        """Generate simulated Dataset

        Args:
            dset:   A dataset containing the data.
        """
        # Initialize variables
        time = list()
        used_satellites = list()
        used_systems = list()
        dset_satellites = list()
        dset_systems = list()
        dset.vars.update(config.files.vars)

        # Select satellites for further orbit comparison
        for satellite in self.satellites:
            if satellite.startswith(tuple(self.systems)):
                used_satellites.append(satellite)
                used_systems.append(satellite[0])  # Needed for generation of Dataset

        # Prepare generation of Dataset
        hour = 0
        minute = 0
        second = 0
        start_time = datetime(
            dset.analysis["rundate"].year,
            dset.analysis["rundate"].month,
            dset.analysis["rundate"].day,
            hour,
            minute,
            second,
        )
        time_to_save = start_time
        num_satellite = len(used_satellites)
        while time_to_save < start_time + timedelta(days=1):
            time.extend([time_to_save] * num_satellite)
            time_to_save += timedelta(seconds=self.sampling_rate)
            dset_satellites.extend(used_satellites)  # All satellites has to be given in each epoch
            dset_systems.extend(used_systems)  # All systems identifiers has to be given in each epoch

        # Generate Dataset
        dset.num_obs = len(time)
        dset.add_time("time", val=time, scale="gps", fmt="datetime", write_level="operational")
        dset.add_text("satellite", val=dset_satellites, write_level="operational")
        dset.add_text("system", val=dset_systems, write_level="operational")
        dset.add_position(
            "site_pos", time=dset.time, val=np.repeat(self.site_pos[None, :], dset.num_obs, axis=0), system="trs"
        )

        # Write Dataset to file
        dset.write_as(stage="setup")
