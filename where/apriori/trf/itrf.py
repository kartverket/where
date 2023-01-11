"""A parser for reading data from ITRF files in SNX format

Description:
------------

Reads station positions and velocities from ITRF files in SNX format. The velocity model is a simple linear offset
based on the reference epoch. For ITRF2014 a model for post seismic deformations are also included. see
:cite:'itrf2014'

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
from where.apriori.trf import TrfFactory
from where.data.position import Position, PositionDelta
from where.data.time import Time
from where.lib import config
from where.lib import exceptions
from where.lib import log
from where import parsers


@plugins.register
class Itrf(TrfFactory):
    """A parser for reading data from ITRF files in SNX format
    """

    file_key_pattern = "trf-itrf_{}"

    def __init__(self, time, version=None):
        """Set up a new Terrestrial Reference Frame object based on the ITRF Sinex files

        For ITRFs from 2014 and later includes post-seismic deformations.

        Args:
            time (Time):    Time epochs for which to calculate positions.
            version (Str):  Version string, can be used to differentiate for instance ITRF2008 from ITRF2014.
        """
        super().__init__(time, version)

        # Parse solution and format from version (solution_format)
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
                log.fatal(f"No ITRF reference frame files found ({file_path}). Find and download one at {url}.")
        self.version = f"{self.solution}_{self.format}"

    @property
    def file_paths(self):
        """File paths used to read sinex data

        TODO: Make this configurable with a files-argument when initializing the factory
        """
        file_vars = dict(version=self.solution)
        paths = {
            self.format: config.files.path(
                self.file_key_pattern.format(self.format), file_vars=file_vars, download_missing=True
            )
        }
        if self.format == "snx":
            if self.solution >= "2014":
                paths.update(
                    dict(
                        soln=config.files.path("trf-itrf_snx_soln", file_vars=file_vars, download_missing=True),
                        psd=config.files.path("trf-itrf_snx_psd", file_vars=file_vars, download_missing=True),
                    )
                )
            else:
                paths.update(
                    dict(soln=config.files.path("trf-itrf_snx_soln", file_vars=file_vars, download_missing=True))
                )
        return paths

    def _read_data(self):
        """Read data needed by this Reference Frame for calculating positions of sites

        Delegates to _read_data_<self.format> to read the actual data.

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        try:
            return getattr(self, f"_read_data_{self.format}")()
        except AttributeError:
            log.fatal(f"Format {self.format!r} is unknown for reference frame {self}")

    def _read_data_ssc(self):
        """Read data needed by ITRF SSC for calculating positions of sites

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        data_trf = parsers.parse_file("trf_ssc", file_path=self.file_paths["ssc"]).as_dict()
        data_ssc = {k.lower(): v for k, v in data_trf.items()}
        return data_ssc

    def _read_data_snx(self):
        """Read data needed by ITRF Sinex for calculating positions of sites

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        paths = self.file_paths
        data_trf = parsers.parse_file("trf_snx", file_path=paths["snx"]).as_dict()
        data_snx = {k.lower(): v for k, v in data_trf.items()}

        for site_key, site_dict in data_snx.items():
            min_soln = min(site_dict["pos_vel"].keys())
            max_soln = max(site_dict["pos_vel"].keys())
            # Change start and end time for first and last solution to allow extrapolation
            site_dict["pos_vel"][min_soln]["start"] = datetime.min
            site_dict["pos_vel"][max_soln]["end"] = datetime.max

        # Post-seismic deformations (for 2014 and later)
        if self.solution >= "2014":
            data_psd = parsers.parse_file("trf_snx_psd", file_path=paths["psd"]).as_dict()
            for site_key, site_dict in data_psd.items():
                site_id = site_key.lower()
                if site_id in data_snx:
                    data_snx[site_id]["psd"] = site_dict["psd"]

        return data_snx

    def _calculate_pos_trs(self, site):
        """Calculate positions for the given time epochs

        The positions are calculated as simple linear offsets based on the reference epoch. If there is a post-seismic
        deformations model for a station the motion due to that model is added to the linear velocity model. Makes sure
        to pick out the correct time interval to use.

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
        pos_trs = Position(np.squeeze(pos), system="trs", ellipsoid=ell, time=self.time)

        # Post-seismic deformations, see Appendix C in :cite:'itrf2014'
        if "psd" in station_info:
            psd = station_info["psd"]
            denu = dict(H=np.zeros(self.time.size), E=np.zeros(self.time.size), N=np.zeros(self.time.size))
            for param in psd.values():
                t_0 = Time(param["epoch"], fmt="datetime", scale="utc")
                delta_t = (self.time - t_0).jd * Unit.day2julian_years
                if isinstance(delta_t, float):
                    delta_t = np.array([delta_t])
                idx = delta_t > 0
                for L in "ENH":
                    aexp = np.array(param.get("AEXP_" + L, list()))
                    texp = np.array(param.get("TEXP_" + L, list()))
                    for a, t in zip(aexp, texp):
                        denu[L][idx] += a * (1 - np.exp(-delta_t[idx] / t))
                    alog = np.array(param.get("ALOG_" + L, list()))
                    tlog = np.array(param.get("TLOG_" + L, list()))
                    for a, t in zip(alog, tlog):
                        denu[L][idx] += a * np.log(1 + delta_t[idx] / t)

            denu = np.vstack((denu["E"], denu["N"], denu["H"])).T

            pos_delta = PositionDelta(np.squeeze(denu), system="enu", ellipsoid=ell, ref_pos=pos_trs, time=self.time)
            pos_trs += pos_delta.trs

        return np.squeeze(pos_trs)
