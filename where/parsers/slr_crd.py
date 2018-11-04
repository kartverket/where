"""A parser for reading SLR data from CRD files

Description:
------------

Reads data from files in the CRD file format as defined in http://ilrs.gsfc.nasa.gov/docs/2009/crd_v1.01.pdf (revision
date October 27, 2009). General information about the format at
http://ilrs.gsfc.nasa.gov/data_and_products/formats/crd.html with links to current specification and data.

"""

# Standard library imports
from datetime import date, datetime, timedelta
import itertools

# External library imports
import numpy as np
from scipy import interpolate

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where.lib import config
from where.lib import files
from where.lib import log
from where.parsers import parser
from where.lib.time import Time
from where.lib.unit import unit


@plugins.register
class SlrCrdParser(parser.Parser):
    """A parser for reading SLR data from CRD files
    """

    def __init__(self, rundate, sat_name):
        super().__init__(rundate)
        self.file_key = "slr_obs_crd"
        self.vars.update(apriori.get_satellite_vars(sat_name))
        self.arc_length = config.tech.arc_length.int

    def read_data(self):
        """Read the data from three monthly datafiles
        """
        files_read = []
        date_to_read = self.rundate - timedelta(days=7)
        while date_to_read < self.rundate + timedelta(days=self.arc_length + 8):
            self.vars.update(config.date_vars(date_to_read))
            file_path = files.path(self.file_key, file_vars=self.vars)
            if file_path not in files_read:
                files_read.append(file_path)
                self.dependencies.append(file_path)
                with files.open(self.file_key, file_vars=self.vars, mode="rt", encoding="latin_1") as fid:
                    self.parse_file(fid)
            date_to_read += timedelta(days=1)

    def setup_parsers(self):
        # Data are organized in sessions
        session_parser = parser.define_parser(
            end_marker=lambda line, _ln, _n: line[0:2].upper() == "H8",
            label=lambda line, _ln: line[0:2].upper(),
            parser_def={
                "H2": {
                    "parser": self.parse_station,
                    "fields": {
                        "cdp": (14, 18),
                        "system_number": (19, 21),
                        "occupancy_number": (22, 24),
                        "epoch_time_scale": (25, 27),
                    },
                },
                "H3": {
                    "parser": self.parse_satellite,
                    "fields": {
                        "satellite_id": (14, 22),
                        "sic": (23, 27),
                        "norad_id": (28, 36),
                        "spacecraft_epoch_time_scale": (37, 38),
                        "target_type": (39, 40),
                    },
                },
                "H4": {
                    "parser": self.parse_session,
                    "fields": {
                        "data_type": (3, 5),
                        "starting_year": (6, 10),
                        "starting_month": (11, 13),
                        "starting_day": (14, 16),
                        "starting_hour": (17, 19),
                        "starting_minute": (20, 22),
                        "starting_second": (23, 25),
                        "ending_year": (26, 30),
                        "ending_month": (31, 33),
                        "ending_day": (34, 36),
                        "ending_hour": (37, 39),
                        "ending_minute": (40, 42),
                        "ending_second": (43, 45),
                        "data_release": (46, 48),
                        "trop_refr_corr_applied": (49, 50),
                        "center_of_mass_corr_applied": (51, 52),
                        "receive_ampl_corr_applied": (53, 54),
                        "station_syst_del_appl": (55, 56),
                        "spacecraft_syst_del_appl": (57, 58),
                        "range_type_indicator": (59, 60),
                        "data_quality_alert": (61, 62),
                    },
                },
                "C0": {"parser": self.parse_config, "fields": ["record_type", "detail_type", "wavelength", None]},
                "C2": {
                    "parser": self.parse_config2,
                    "fields": [
                        "record_type",
                        "detail_type",
                        "detector_configuration_id",
                        "detector_type",
                        "applicable_wavelength",
                        "quantum_efficiency",
                        "applied_voltage",
                        "dark_count",
                        "output_pulse_type",
                        "output_pulse_width",
                        "spectral_filter",
                        None,
                    ],
                },
                "11": {"parser": self.parse_obs_normal, "fields": [None, "second", "time_of_flight"]},
                "20": {
                    "parser": self.parse_meteorology,
                    "fields": [None, "seconds", "pressure", "temperature", "humidity"],
                },
            },
        )
        return itertools.repeat(session_parser)

    def parse_station(self, line, cache):
        cache["station"] = line.pop("cdp")
        if "station_{}".format(cache["station"]) in self.data:
            return

        self.data["station_{}".format(cache["station"])] = line

    def parse_satellite(self, line, cache):
        cache["satellite"] = line.pop("satellite_id")
        self.data["satellite_{}".format(cache["satellite"])] = line

    def parse_session(self, line, cache):
        cache["obs_date_first"] = datetime(
            int(line["starting_year"]),
            int(line["starting_month"]),
            int(line["starting_day"]),
            int(line["starting_hour"]),
            int(line["starting_minute"]),
            int(line["starting_second"]),
        )
        cache["starting_hour"] = int(line["starting_hour"])
        cache["obs_date"] = date(int(line["starting_year"]), int(line["starting_month"]), int(line["starting_day"]))
        cache["session"] = "{station}/{satellite}/{obs_date_first}".format(**cache)
        cache["data_type"] = line.pop("data_type")

    def parse_config(self, line, cache):
        cache["wavelength"] = line.pop("wavelength")

    def parse_config2(self, line, cache):
        cache["detector_type"] = line.pop("detector_type") if "detector_type" in line else ""

    def parse_obs_normal(self, line, cache):
        """Parse the observation line, normal point data
        """
        obs_time = (cache["obs_date"] - self.rundate).total_seconds() + float(line["second"])
        if float(line["second"]) < 43000.0 and cache["starting_hour"] > 12:
            # Suspect day-shift in session:
            obs_time = obs_time + timedelta(days=1).total_seconds()

        # Use one week (or arc_length) of data starting with rundate:
        if 0 <= obs_time < timedelta(days=self.arc_length).total_seconds():
            # Meta information about observation
            obs_meta = {
                "time": obs_time * unit.sec2day,
                "station": cache["station"],
                "satellite": cache["satellite"],
                "session": cache["session"],
            }
            for field, value in obs_meta.items():
                self.meta.setdefault(field, list()).append(value)

            # Float data concerning the observation
            obs = {
                "obs_time": obs_time,
                "time_of_flight": float(line["time_of_flight"]),
                "wavelength": float(cache["wavelength"]),
                "data_type": float(cache["data_type"]),
            }
            for field, value in obs.items():
                self.data.setdefault("obs", dict()).setdefault(field, list()).append(value)

            # Text data concerning the observation
            obs_str = {"detector_type": cache.get("detector_type", "")}
            for field, value in obs_str.items():
                self.data.setdefault("obs_str", dict()).setdefault(field, list()).append(value)

    def parse_meteorology(self, line, cache):
        """Parse the meteorological record
        """
        met_secs = (cache["obs_date"] - self.rundate).total_seconds() + float(line["seconds"])
        if float(line["seconds"]) < 43000 and cache["starting_hour"] > 12:  # Suspect day-shift in session:
            met_secs += 86400

        met_data = {
            "time": met_secs * unit.sec2day,
            "temperature": float(line["temperature"]),
            "humidity": float(line["humidity"]),
            "pressure": float(line["pressure"]),
        }
        for field, value in met_data.items():
            self.data.setdefault("met_" + cache["station"], dict()).setdefault(field, list()).append(value)

        met_line = (float(line["temperature"]), float(line["humidity"]), float(line["pressure"]))
        self.data.setdefault("met", {}).setdefault(cache["station"], {})
        self.data["met"][cache["station"]][met_secs] = met_line

    #
    # CALCULATORS for calculating other necessary data
    #
    def setup_calculators(self):
        return [self.interpolate_meteorological_data]

    def interpolate_meteorological_data(self):
        """Calculate temperature, humidity and pressure at observation epochs

        Meteorological data are calculated at observation epochs by interpolating in the data given on the observation
        file for each station.

        Missing meteorological data are currently not handled.
        """
        obs_time = self.meta["time"]
        max_obs = max(obs_time)
        min_obs = min(obs_time)

        for field, station in [(f, f[4:]) for f in self.data.keys() if f.startswith("met_")]:
            log.debug("Meteorological data available for station {}", station)

            met_time = self.data[field].pop("time")
            met_time_sorted = sorted(met_time)
            max_met_time = met_time_sorted[-1]
            min_met_time = met_time_sorted[0]
            max_idx = met_time.index(max_met_time)
            min_idx = met_time.index(min_met_time)
            for met_type in self.data[field].keys():
                temp_array = np.zeros(len(obs_time))
                for i in range(0, len(obs_time)):
                    # Extrapolating one hour before and after available met data, other missing met data is set to zero
                    if min_met_time <= obs_time[i] <= max_met_time:
                        temp_array[i] = interpolate.interp1d(met_time, self.data[field][met_type])(obs_time[i])
                    elif min_met_time - 1 / 24 < obs_time[i] < min_met_time:
                        temp_array[i] = self.data[field][met_type][min_idx]
                    elif max_met_time < obs_time[i] < max_met_time + 1 / 24:
                        temp_array[i] = self.data[field][met_type][max_idx]
                self.data[field][met_type] = temp_array

    def write_to_dataset(self, dset):
        """Store SLR data in a dataset

        Args:
           dset_out: The Dataset where data are stored.
        """
        dset.num_obs = len(self.meta["time"])
        dset.add_time(
            "time", val=Time(val=self.rundate.isoformat()).mjd, val2=self.meta.pop("time"), scale="utc", format="mjd"
        )
        for field, value in self.meta.items():
            dset.add_text(field, val=value)

        # Positions
        trf = apriori.get("trf", time=dset.time)
        for station in dset.unique("station"):
            trf_site = trf[station]
            station_pos = trf_site.pos.itrs
            log.debug(
                "Station position for {} ({}) according to ITRF is (x,y,z) = {}",
                station,
                trf_site.name,
                station_pos.mean(axis=0),
            )
            domes = trf_site.meta["domes"]

            if False:  # TODO: Add these missing stations to trf-file
                domes = "00000"
                log.warn("No information about station {} on ITRF file", station)
                if station == "7407":
                    station_pos = np.repeat([[4119502.13, -4553595.23, -1722855.13]], dset.num_obs, axis=0)
                elif station == "1889":
                    station_pos = np.repeat([[3451136.221, 3060335.064, 4391970.241]], dset.num_obs, axis=0)
                elif station == "1888":
                    station_pos = np.repeat([[2730139.097, 1562328.629, 5529998.585]], dset.num_obs, axis=0)
                elif station == "1891":
                    station_pos = np.repeat([[-968340.32, 3794415.10, 5018178.10]], dset.num_obs, axis=0)
                elif station == "1887":
                    station_pos = np.repeat([[2001873.3346, 3987633.3547, 4542477.6716]], dset.num_obs, axis=0)
                elif station == "1886":
                    station_pos = np.repeat([[3466773.394, 3059757.864, 4381456.782]], dset.num_obs, axis=0)
                elif station == "1874":
                    station_pos = np.repeat([[2844591.641, 2161111.997, 5266356.839]], dset.num_obs, axis=0)
                elif station == "1890":
                    station_pos = np.repeat([[-838299.699, 3865738.865, 4987640.921]], dset.num_obs, axis=0)
                else:
                    log.error("Unknown station {}", station)
                    station_pos = np.zeros((dset.num_obs, 3))
                log.warn("Using coordinates {} for {}", np.mean(station_pos, axis=0), station)

            self.data["pos_" + station] = station_pos
            self.data["station-other_" + station] = dict(domes=domes, cdp=station, site_id=station)

        dset.add_position(
            "site_pos", time="time", itrs=np.array([self.data["pos_" + s][idx] for idx, s in enumerate(dset.station)])
        )

        # Station data
        sta_fields = set().union(*[v.keys() for k, v in self.data.items() if k.startswith("station_")])
        for field in sta_fields:
            dset.add_float(field, val=np.array([float(self.data["station_" + s][field]) for s in dset.station]))
        sta_fields = set().union(*[v.keys() for k, v in self.data.items() if k.startswith("station-other_")])
        for field in sta_fields:
            dset.add_text(field, val=[self.data["station-other_" + s][field] for s in dset.station])

        # Satellite data
        sat_fields = set().union(*[v.keys() for k, v in self.data.items() if k.startswith("satellite_")])
        for field in sat_fields:
            dset.add_float(field, val=np.array([float(self.data["satellite_" + s][field]) for s in dset.satellite]))

        # Observations
        for field, values in self.data["obs"].items():
            dset.add_float(field, val=np.array(values))

        for field, values in self.data["obs_str"].items():
            dset.add_text(field, val=values)

        # Meterological data
        met_fields = set().union(*[v.keys() for k, v in self.data.items() if k.startswith("met_")])
        for field in met_fields:
            dset.add_float(field, val=np.diag([self.data["met_" + s][field] for s in dset.station]))
