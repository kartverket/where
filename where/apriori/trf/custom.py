"""A TRF Factory for creating sites with positions from a custom config file

Description:
------------

Reads station positions and velocities from a custom config file. No velocities are available in the file, so the same
position is used for all time epochs.

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.config.config import Configuration
from midgard.dev import plugins
from midgard.files import dependencies
from midgard.math import ellipsoid
from midgard.math.unit import Unit

# Where imports
from where.data.position import Position
from where.data.time import Time
from where.apriori.trf import TrfFactory
from where.lib import config


@plugins.register
class CustomTrf(TrfFactory):
    """A factory for using positions from the custom TRF-config file.
    """

    def _read_data(self):
        """Read data needed by this Reference Frame for calculating positions of sites

        Delegates to _read_data_<self.version> to read the actual data.

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        trf = Configuration("trf_custom")
        trf.profiles = config.analysis.get("pipeline", value=self.version, default="").list
        trf_path = config.files.path("trf-custom")
        trf_local_path = config.files.path("trf-custom_local")

        trf.update_from_file(trf_path)
        dependencies.add(trf_path, label="trf")
        if trf_local_path.exists():
            trf.update_from_file(trf_local_path)
            dependencies.add(trf_local_path, label="trf")

        data = dict()
        for section in trf.sections:
            info = {k: v for k, v in section.as_dict().items() if not k == "pos_itrs" and not k == "vel_itrs"}
            info["pos"] = np.array(section.pos_itrs.list, dtype=float)
            if "vel_itrs" in section and "ref_epoch" in section:
                info["vel"] = np.array(section.vel_itrs.list, dtype=float)
            data[section.name] = info

        return data

    def _calculate_pos_trs(self, site):
        """Calculate positions for the given time epochs

        There are no velocities available, so same position is returned for all time epochs

        Args:
            site (String):  Key specifying which site to calculate position for, must be key in self.data.

        Returns:
            Array:  Positions, one 3-vector for each time epoch.
        """
        station_info = self.data[site]
        ell = ellipsoid.get(config.tech.reference_ellipsoid.str.upper())
        if "vel" in station_info and "ref_epoch" in station_info:
            ref_epoch = Time(float(info["ref_epoch"]), scale="utc", fmt="decimalyear")
            ref_pos = station_info.pop("pos")
            ref_vel = station_info.pop("vel")
            interval_years = (self.time - ref_epoch).jd * Unit.day2julian_years
            if isinstance(interval_years, float):
                interval_years = np.array([interval_years])

            pos = ref_pos + interval_years[:, None] * ref_vel[None, :]
        else:
            pos = station_info.pop("pos")[None, :].repeat(self.time.size, axis=0)

        pos_trs = Position(pos, system="trs", ellipsoid=ell, time=self.time)
        return np.squeeze(pos_trs)
