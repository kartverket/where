"""Read and parse the Rinex Doris 3.0 format.

Description:
------------

The Rinex Doris 3.0 format is specified in [1]. This is a simple extension of the Rinex (The Receiver Independent
Exchange) 3.0 format [2], originally developed by the Astronomical Institute of the University of Berne for the easy
exchange of GPS data [2, page 4].

References:
-----------

.. [1] E. Lourme, "Rinex Doris 3.0". Dated 11/05/2010.
   Available at ftp://ftp.ids-doris.org/pub/ids/data/RINEX_DORIS.pdf

.. [2] W. Gurtner, and Estey. L., "RINEX. The Receiver Independent
   Exchange Format. Version 3.00". Dated November 28, 2007.
   Available at ftp://igs.org/pub/data/format/rinex300.pdf

"""

# Standard library imports
import itertools

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where.lib import config
from where.lib import log
from where.parsers import parser


@plugins.register
class DorisRinex3Parser(parser.Parser):
    """A parser for reading DORIS data from Rinex 3 files
    """

    def __init__(self, rundate, satellite):
        super().__init__(rundate)
        self.file_key = "doris_obs_rinex3"
        self.vars["sat_name"] = satellite
        self.vars["sat_shortname"] = "ja2"  # TODO

        # Read arc-length from config file
        self.arc_length = config.tech.arc_length.int

        #    def read_data(self):
        #        """Read input data from all session files within the arc length  TODO
        #        """

    def setup_parsers(self):
        header_parser = parser.define_parser(
            end_marker=lambda line, _ln, _n: line[60:73] == "END OF HEADER",
            label=lambda line, _ln: line[60:].strip(),
            parser_def={
                "REC # / TYPE / VERS": {"parser": self._parse_default, "fields": {"receiver": (20, 40)}},
                "ANT # / TYPE": {"parser": self._parse_default, "fields": {"antenna_type": (20, 40)}},
                "APPROX POSITION XYZ": {
                    "parser": self._parse_xyz,
                    "fields": {"antenna_pos_x": (0, 14), "antenna_pos_y": (14, 28), "antenna_pos_z": (28, 42)},
                },
                "CENTER OF MASS: XYZ": {
                    "parser": self._parse_xyz,
                    "fields": {
                        "center_of_mass_x": (0, 14),
                        "center_of_mass_y": (14, 28),
                        "center_of_mass_z": (28, 42),
                    },
                },
                "SYS / SCALE FACTOR": {
                    "parser": self._parse_default,
                    "fields": {"scale_factor": (2, 6), "fields": (10, 58)},
                },
                "L2 / L1 DATE OFFSET": {"parser": self._parse_float, "fields": {"date_offset": (3, 17)}},
                "STATION REFERENCE": {
                    "parser": self._parse_station,
                    "fields": {
                        "station_key": (0, 3),
                        "site_id": (5, 9),
                        "site_name": (10, 40),
                        "domes": (40, 50),
                        "type": (51, 52),
                        "freq_shift": (53, 56),
                    },
                },
                "TIME REF STATION": {
                    "parser": self._parse_timeref_station,
                    "fields": {"station_key": (0, 3), "bias": (5, 19), "shift": (21, 35)},
                },
                "TIME REF STAT DATE": {
                    "parser": self._parse_float,
                    "fields": {
                        "timeref_year": (0, 6),
                        "timeref_month": (6, 12),
                        "timeref_day": (12, 18),
                        "timeref_hour": (18, 24),
                        "timeref_minute": (24, 30),
                        "timeref_second": (30, 43),
                    },
                },
            },
        )

        """Definition of the observation data of the input file.

        This dictionary defines how to parse the data part of the input file, and is based on the "Rinex Doris
        3.0"-document (Appendix A).
        """
        obs_parser = parser.define_parser(
            end_marker=lambda _l, _ln, next_line: next_line.startswith(">"),
            label=lambda line, _ln: line[0:1],
            parser_def={
                ">": {
                    "parser": self._parse_epoch,
                    "fields": {
                        "year": (2, 6),
                        "month": (7, 9),
                        "day": (10, 12),
                        "hour": (13, 15),
                        "minute": (16, 18),
                        "second": (19, 21),
                        "fracsec": (21, 31),
                        "epoch_flag": (33, 34),
                        "numrecords": (34, 37),
                        "clock_offset": (43, 56),
                        "clock_offset_flag": (57, 58),
                    },
                },
                "D": {
                    "parser": self._parse_observation_line_1,
                    "fields": {
                        "station_key": (0, 3),
                        "l1": (3, 17),
                        "l1_flag": (17, 19),
                        "l2": (19, 33),
                        "l2_flag": (33, 35),
                        "c1": (35, 49),
                        "c1a_flag": (49, 50),
                        "c1b_flag": (50, 51),
                        "c2": (51, 65),
                        "c2a_flag": (65, 66),
                        "c2b_flag": (66, 67),
                        "w1": (67, 81),
                        "w1_flag": (81, 83),
                    },
                },
                " ": {
                    "parser": self._parse_observation_line_2,
                    "fields": {
                        "w2": (3, 17),
                        "w2_flag": (17, 19),
                        "f": (19, 33),
                        "f_flag": (33, 35),
                        "p": (35, 49),
                        "p_flag": (49, 51),
                        "t": (51, 65),
                        "t_flag": (65, 67),
                        "h": (67, 81),
                        "h_flag": (81, 83),
                    },
                },
            },
        )
        return itertools.chain([header_parser], itertools.repeat(obs_parser))

    def _parse_default(self, line, _):
        """A general default parser: Adds all fields in `line` to `data`.

        Parameters
        ----------
        data : dict
            A dictionary containing all data already parsed, this dictionary is being updated by this parse function.
        line : dict
            A dictionary containing the data to be parsed.
        _ : None, optional
            Ignored. Needed because when parsing observation data, some cached data is passed along.
        """
        self.data.update(line)

    def _parse_float(self, line, _):
        """A general parser for floating point numbers.

        See Also
        ----------
        _parse_default
        """
        self._parse_default({k: float(line[k]) for k in line if line[k]}, _)

    def _parse_xyz(self, line, _):
        """A general parser for xyz positions.

        See Also
        ----------
        _parse_default
        """
        xyz_fields = {f[:-2] for f in line.keys() if f[-2] == "_" and f[-1] in "xyz"}
        for field in xyz_fields:
            x, y, z = [float(line.pop(field + "_" + xyz)) for xyz in "xyz"]
            self.data[field] = np.array([x, y, z])
        self._parse_default(line, _)

    def _parse_station(self, line, _):
        """Parse information about a given station.

        Checks the station information against the known data, and warns about differences.

        See Also
        --------
        _parse_default

        """
        # Antenna of station is identified by the last character of the station code
        antenna_map = {"A": 0, "B": 1, "C": 2}  # TODO: Use antenna names

        # Create a station dict to store information about the station
        station_key = line.pop("station_key")
        site_id = line["site_id"]
        line["site_name"] = line["site_name"].title()
        station = dict(**line)

        # Add extra information to station dictionary
        station["antenna_flag"] = antenna_map[site_id[-1]]
        station["timeref_bias"] = 0.0  # Default value, possibly overwritten in _parse_timeref_station
        station["timeref_shift"] = 0.0  # Default value, possibly overwritten in _parse_timeref_station

        # Add all information (including bookkeeping) to meta dictionary
        self.meta.setdefault("station_map", dict())[station_key] = site_id
        self.meta.setdefault("station_info", dict())[site_id] = station
        self.meta.setdefault("station_list", set()).add(site_id)

    def _parse_timeref_station(self, line, _):
        """Parse information about the reference stations.

        Adds the information to the relevant station

        See Also
        --------
        _parse_default
        """
        station_key = line.pop("station_key")
        site_id = self.meta["station_map"][station_key]
        # self._parse_float(line, line)

        self.meta["station_info"][site_id]["timeref_bias"] = line["bias"]
        self.meta["station_info"][site_id]["timeref_shift"] = line["shift"]

    def _parse_epoch(self, line, cache):
        """Parse the epoch of an observation.

        See Also
        --------
        _parse_default
        """
        cache["obs"] = {
            "time": "{year}-{month}-{day}T{hour}:{minute}:{second}".format(**line),
            "fracsec": line["fracsec"],
            "clock_offset": float(line["clock_offset"]),
        }

    def _parse_observation_line_1(self, line, cache):
        """Parse the first line of an observation for a station.

        See Also
        --------
        _parse_default
        """
        station_key = line.pop("station_key")
        cache["obs"]["station"] = self.meta["station_map"][station_key]
        cache["obs"].update({k: float(v) for k, v in line.items() if not k.endswith("_flag")})

    def _parse_observation_line_2(self, line, cache):
        """Parse the second line of an observation for a station.

        See Also
        --------
        _parse_default
        """
        cache["obs"].update({k: float(v) for k, v in line.items() if not k.endswith("_flag")})
        for field, value in cache["obs"].items():
            self.data.setdefault("obs", dict()).setdefault(field, list()).append(value)

    #
    # WRITE DATA
    #
    def write_to_dataset(self, dset):
        """Store DORIS data in a dataset

        Args:
            dset: The Dataset where data are stored.
        """
        dset.num_obs = len(self.data["obs"]["time"])
        dset.add_time("time", val=self.data["obs"].pop("time"), scale="utc", format="isot")
        dset.add_text("station", val=self.data["obs"].pop("station"))

        for field, value in self.data["obs"].items():
            dset.add_float(field, val=np.array(value))

        # Station data
        sta_fields = set().union(*[v.keys() for v in self.meta["station_info"].values()])
        for field in sta_fields:
            dset.add_text(field, val=[self.meta["station_info"][s][field] for s in dset.station])

        # Station positions
        site_pos = np.zeros((dset.num_obs, 3))
        trf = apriori.get("trf", time=dset.time)
        for site in dset.unique("station"):
            idx = dset.filter(station=site)
            site_pos[idx, :] = trf[site].pos.trs[idx, :]
            log.debug(f"Using position {np.mean(site_pos[idx, :], axis=0)} for {site!r}")
        dset.add_position("site_pos", time="time", itrs=site_pos)

        # Satellite
        dset.add_text("satellite", val=[self.vars["sat_name"]] * dset.num_obs)
