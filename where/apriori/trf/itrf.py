"""A parser for reading data from ITRF files in SNX format

Description:
------------

Reads station positions and velocities from ITRF files in SNX format. The velocity model is a simple linear offset
based on the reference epoch. For ITRF2014 a model for post seismic deformations are also included. see
:cite:'itrf2014'

References:
-----------




"""

# External library imports
import numpy as np

# Where imports
from where.apriori import trf
from where.lib import files
from where.lib import log
from where.lib import plugins
from where import parsers
from where.lib import rotation
from where.ext import sofa_wrapper as sofa
from where.lib.time import Time
from where.lib.unit import unit


@plugins.register
class Itrf(trf.TrfFactory):
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
            candidates = files.glob_variable(self.file_key_pattern.format(self.format), "version", r"[^._]*")
            try:
                self.solution = max(candidates)
            except ValueError:
                log.fatal("No itrf reference frame files found")
        self.version = "{}_{}".format(self.solution, self.format)

    @property
    def file_paths(self):
        """File paths used to read sinex data

        TODO: Make this configurable with a files-argument when initializing the factory
        """
        file_vars = dict(version=self.solution)
        paths = {
            self.format: files.path(
                self.file_key_pattern.format(self.format), file_vars=file_vars, download_missing=True
            )
        }
        if self.format == "snx":
            if self.solution >= "2014":
                paths.update(
                    dict(
                        soln=files.path("trf-itrf_snx_soln", file_vars=file_vars, download_missing=True),
                        psd=files.path("trf-itrf_snx_psd", file_vars=file_vars, download_missing=True),
                    )
                )
            else:
                paths.update(dict(soln=files.path("trf-itrf_snx_soln", file_vars=file_vars, download_missing=True)))
        return paths

    def _read_data(self):
        """Read data needed by this Reference Frame for calculating positions of sites

        Delegates to _read_data_<self.format> to read the actual data.

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        try:
            return getattr(self, "_read_data_{}".format(self.format))()
        except AttributeError:
            log.fatal("Format '{}' is unknown for reference frame {}", self.format, self)

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
        paths = self.file_paths
        data_trf = parsers.parse_file("trf_snx", file_path=paths["snx"]).as_dict()
        data_snx = {k.lower(): v for k, v in data_trf.items()}
        # Time epoch intervals are in a separate file
        data_soln = parsers.parse_file("trf_snx_soln", file_path=paths["soln"]).as_dict()
        for site_key, site_dict in data_soln.items():
            site_id = site_key.lower()
            for soln, interval in site_dict.items():
                if soln in data_snx[site_id]["pos_vel"]:
                    data_snx[site_id]["pos_vel"][soln]["start"] = interval["start"]
                    data_snx[site_id]["pos_vel"][soln]["end"] = interval["end"]
                elif soln - 1 in data_snx[site_id]["pos_vel"]:
                    # copy previous solution for extrapolation
                    data_snx[site_id]["pos_vel"][soln] = data_snx[site_id]["pos_vel"][soln - 1].copy()
                    data_snx[site_id]["pos_vel"][soln]["start"] = interval["start"]
                    data_snx[site_id]["pos_vel"][soln]["end"] = interval["end"]

        # Post-seismic deformations (for 2014 and later)
        if self.solution >= "2014":
            data_psd = parsers.parse_file("trf_snx_psd", file_path=paths["psd"]).as_dict()
            for site_key, site_dict in data_psd.items():
                site_id = site_key.lower()
                if site_id in data_snx:
                    data_snx[site_id]["psd"] = site_dict["psd"]

        return data_snx

    def _calculate_pos_itrs(self, site):
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
        ref_epoch = Time(station_info["ref_epoch"], scale="utc")

        pos = np.zeros((self.time.size, 3))
        for pv in station_info["pos_vel"].values():
            idx = np.logical_and(self.time.utc.datetime >= pv["start"], self.time.utc.datetime < pv["end"])
            if self.time.size == 1:
                idx = np.array([idx])
            if not any(idx):
                continue
            ref_pos = np.array([pv["STAX"], pv["STAY"], pv["STAZ"]])
            ref_vel = np.array([pv["VELX"], pv["VELY"], pv["VELZ"]])
            interval_years = (self.time - ref_epoch).jd * unit.day2julian_years
            if isinstance(interval_years, float):
                interval_years = np.array([interval_years])
            pos[idx, :] = ref_pos + interval_years[idx, None] * ref_vel[None, :]

        # Post-seismic deformations, see Appendix C in :cite:'itrf2014'
        if "psd" in station_info:
            llh = sofa.vectorized_llh(pos)
            psd = station_info["psd"]
            denu = dict(H=np.zeros(self.time.size), E=np.zeros(self.time.size), N=np.zeros(self.time.size))
            for param in psd.values():
                t_0 = Time(param["epoch"], format="datetime", scale="utc")
                delta_t = (self.time - t_0).jd * unit.day2julian_years
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

            rot = rotation.enu2trf(llh[:, 0], llh[:, 1])
            denu = np.vstack((denu["E"], denu["N"], denu["H"])).T
            dxyz = (rot @ denu[:, :, None])[:, :, 0]
            pos += dxyz

        if self.time.size == 1:
            pos = pos[0, :]
        return pos
