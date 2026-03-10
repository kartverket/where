"""A parser for reading data from VTRF files in SNX format

Description:
------------

Reads station positions and velocities from VTRF files in SNX format. The velocity model is a simple linear offset
based on the reference epoch. VTRF is normally released four times per year by the IVS combination center.

References:
-----------




"""
from datetime import datetime

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math import ellipsoid
from midgard.math.unit import Unit

# Where imports
from where.data.position import Position
from where.data.time import Time
from where.apriori import trf
from where.lib import config
from where.lib import exceptions
from where.lib import log
from where import parsers


@plugins.register
class Vtrf(trf.TrfFactory):
    """A parser for reading data from ITRF files in SNX format
    """

    file_key_pattern = "trf-vtrf_{}"

    def __init__(self, time, version=None):
        """Set up a new Terrestrial Reference Frame object based on the VTRF Sinex files

        Args:
            time (Time):    Time epochs for which to calculate positions.
            version (Str):  Version string, can be used to differentiate for instance VTRF2015d from VTRF2016b. By
                            adding a _snx or _ssc suffix to the version number the format can be specificed.
        """
        super().__init__(time, version)

        # Parse solution and format from version (year_format)
        version = "last" if version is None else version  # 'last' is default solution
        self.solution, _, fmt = version.partition("_")
        self.format = fmt if fmt else "snx"  # Sinex is default format
        if self.solution == "last":
            file_key = self.file_key_pattern.format(self.format)
            candidates = config.files.glob_variable(file_key, "version", r"[^._]*")
            try:
                self.solution = max(candidates)
            except ValueError:
                url = config.files.url(file_key).base
                file_path = config.files.path(file_key)
                log.fatal(f"No VTRF reference frame files found ({file_path}). Find and download one at {url}.")
        self.version = f"{self.solution}_{self.format}"

    @property
    def file_paths(self):
        """File paths used to read VTRF data (dependent on format)

        TODO: Make this configurable with a files-argument when initializing the factory
        """
        file_vars = dict(version=self.solution) if self.solution else None
        return {
            self.format: config.files.path(
                self.file_key_pattern.format(self.format), file_vars=file_vars, download_missing=True
            )
        }

    def _read_data(self):
        """Read data needed by this Reference Frame for calculating positions of sites

        Delegates to _read_data_<self.format> to read the actual data.

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        try:
            return getattr(self, f"_read_data_{self.format}")()
        except AttributeError:
            log.fatal(f"Format {self.format} is unknown for reference frame {self}")

    def _read_data_ssc(self):
        """Read data needed by ITRF SSC for calculating positions of sites

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        return parsers.parse_file("trf_ssc", file_path=self.file_paths["ssc"]).as_dict()

    def _read_data_snx(self):
        """Read data needed by ITRF Sinex for calculating positions of sites

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        data = parsers.parse_file("trf_snx", file_path=self.file_paths["snx"]).as_dict()
        # A station may be defined in the sinex file but not in included in the solution/estimate block.
        # These stations are discarded
        missing_posvel = [key for key in data.keys() if "ref_epoch" not in data[key]]
        for site_key in missing_posvel:
            del data[site_key]

        for site_key, site_dict in data.items():
            min_soln = min(site_dict["pos_vel"].keys())
            max_soln = max(site_dict["pos_vel"].keys())
            # Change start and end time for first and last solution to allow extrapolation
            site_dict["pos_vel"][min_soln]["start"] = datetime.min
            site_dict["pos_vel"][max_soln]["end"] = datetime.max

        return data

    def _calculate_pos_trs(self, site):
        """Calculate positions for the given time epochs

        The positions are calculated as simple linear offsets based on the reference epoch.

        Args:
            site (String):    Key saying which site to calculate position for.

        Returns:
            Array:  Positions, one 3-vector for each time epoch.
        """
        station_info = self.data[site]
        ref_epoch = Time(station_info["ref_epoch"], scale="utc", fmt="datetime")

        pos = np.zeros((self.time.size, 3))
        for pv in station_info["pos_vel"].values():
            idx = np.logical_and(self.time.utc.datetime >= pv["start"], self.time.utc.datetime < pv["end"])
            if idx.ndim == 0:
                idx = np.array([idx])
            if not any(idx):
                continue
            ref_pos = np.array([pv["STAX"], pv["STAY"], pv["STAZ"]])
            ref_vel = np.array([pv["VELX"], pv["VELY"], pv["VELZ"]])
            interval_years = (self.time - ref_epoch).jd * Unit.day2julian_years
            if isinstance(interval_years, float):
                interval_years = np.array([interval_years])
            pos[idx, :] = ref_pos + interval_years[idx, None] * ref_vel[None, :]

        if not pos.any():
            # All positions are zero
            raise exceptions.MissingDataError(f"Position for {site} is not well defined in {self}")

        ell = ellipsoid.get(config.tech.reference_ellipsoid.str.upper())
        pos_trs = Position(pos, system="trs", ellipsoid=ell, time=self.time)

        return np.squeeze(pos_trs)
