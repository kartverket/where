"""A parser for reading SLR data from CRD files

Description:
------------

Reads data from files in the CRD file format as defined in http://ilrs.gsfc.nasa.gov/docs/2009/crd_v1.01.pdf (revision
date October 27, 2009). General information about the format at
http://ilrs.gsfc.nasa.gov/data_and_products/formats/crd.html with links to current specification and data.

"""

# Standard library imports
from datetime import datetime, timedelta
import itertools

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.parsers._parser_chain import ParserDef, ChainParser


@plugins.register
class SlrCrdParser(ChainParser):
    """A parser for reading SLR data from CRD files
    """

    def setup_parser(self):
        # Data are organized in sessions
        session_parser = ParserDef(
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
                "11": {
                    "parser": self.parse_obs_normal,
                    "fields": [None, "second", "time_of_flight", "epoch_event", "number_of_raw", "bin_rms"],
                },
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
        cache["obs_date"] = datetime(
            int(line["starting_year"]), int(line["starting_month"]), int(line["starting_day"])
        )
        cache["session"] = "{station}/{satellite}/{obs_date_first}".format(**cache)
        cache["data_type"] = line.pop("data_type")

    def parse_config(self, line, cache):
        cache["wavelength"] = line.pop("wavelength")

    def parse_config2(self, line, cache):
        cache["detector_type"] = line.pop("detector_type") if "detector_type" in line else ""

    def parse_obs_normal(self, line, cache):
        """Parse the observation line, normal point data
        """
        obs_sec = float(line["second"])
        obs_time = cache["obs_date"] + timedelta(seconds=obs_sec)
        if obs_sec < 43000.0 and cache["starting_hour"] > 12:
            # Suspect day-shift in session:
            obs_time = obs_time + timedelta(days=1)

        # Meta information about observation
        obs_meta = {
            "obs_time": obs_time,
            "obs_date": cache["obs_date"],
            "obs_sec": obs_sec,
            "station": cache["station"],
            "satellite": cache["satellite"],
            "session": cache["session"],
            "bin_rms": float(line["bin_rms"]),
        }

        for field, value in obs_meta.items():

            self.data.setdefault("meta", dict()).setdefault(field, list()).append(value)

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
        met_secs = float(line["seconds"])
        if float(line["seconds"]) < 43000 and cache["starting_hour"] > 12:  # Suspect day-shift in session:
            met_secs += 86400
        met_time = cache["obs_date"] + timedelta(seconds=met_secs)

        met_data = {
            "met_time": met_time,
            "temperature": float(line["temperature"]),
            "humidity": float(line["humidity"]),
            "pressure": float(line["pressure"]),
        }
        for field, value in met_data.items():
            self.data.setdefault("met_" + cache["station"], dict()).setdefault(field, list()).append(value)
