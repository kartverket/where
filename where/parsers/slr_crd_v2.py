"""A parser for reading SLR data from CRD files

Description:
------------

Reads data from files in the CRD file format version 2 as defined in
https://ilrs.gsfc.nasa.gov/docs/2019/crd_v2.01.pdf (revision date October 19, 2019).

General information about the format at
https://ilrs.gsfc.nasa.gov/data_and_products/formats/crd.html
with links to current specification and data.

"""

# Standard library imports
from datetime import datetime, timedelta
import itertools

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_chain import ParserDef, ChainParser


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
                    "fields": [
                        "cdp",
                        "system_number",
                        "occupancy_number",
                        "epoch_time_scale",
                        "station_network"
                    ],
                },
                "H3": {
                    "parser": self.parse_satellite,
                    "fields": [
                        "satellite_name",
                        "satellite_id",
                        "sic",
                        "norad_id",
                        "spacecraft_epoch_time_scale",
                        "target_type",
                        "target_class",
                        "target_location"
                    ],
                },
                "H4": {
                    "parser": self.parse_session,
                    "fields": [
                        "data_type",
                        "starting_year",
                        "starting_month",
                        "starting_day",
                        "starting_hour",
                        "starting_minute",
                        "starting_second",
                        "ending_year",
                        "ending_month",
                        "ending_day",
                        "ending_hour",
                        "ending_minute",
                        "ending_second",
                        "data_release",
                        "trop_refr_corr_applied",
                        "center_of_mass_corr_applied",
                        "receive_ampl_corr_applied",
                        "station_syst_del_appl",
                        "spacecraft_syst_del_appl",
                        "range_type_indicator",
                        "data_quality_alert"
                    ],
                },
                "H5": {
                    "parser": self.parse_prediction,
                    "fields": [
                        "prediction_type",
                        "date_and_time",
                        "prediction_provider",
                        "sequence_number"
                    ],
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
        cache["satellite"] = line.pop("satellite_id").lstrip("0")
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
        if "data_release" in cache:
            if not cache["data_release"] == 0:
                log.dev("TODO: Replacement release of data, do something")

    def parse_config(self, line, cache):
        cache["wavelength"] = line.pop("wavelength")

    def parse_config2(self, line, cache):
        cache["detector_type"] = line.pop("detector_type") if "detector_type" in line else ""

    def parse_obs_normal(self, line, cache):
        """Parse the observation line, normal point data
        """
        obs_sec = float(line["second"])
        obs_date = cache["obs_date"]

        if obs_sec < 43000.0 and cache["starting_hour"] > 12:
            # Suspect day-shift in session:
            obs_date = obs_date + timedelta(days=1)

        obs_time = obs_date + timedelta(seconds=obs_sec)

        # Meta information about observation
        obs_meta = {
            "obs_time": obs_time,
            "obs_date": obs_date,
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
