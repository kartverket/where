"""A parser for reading SLR center of mass corrections from file

Description:
------------

Asdf

"""

# Standard library imports
from datetime import datetime
import itertools

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import log
from where.parsers import parser
from where.lib.unit import unit


DEFAULT_COM = {"lageos": 0.251, "etalon": 0.65}


@plugins.register
class SlrComParser(parser.ParserDict):
    """A parser for reading SLR center of mass corrections from files
    """

    def __init__(self, sat_name):
        super().__init__()
        self.vars["satellite"] = sat_name.lower().rstrip("12")
        self.file_key = "center_of_mass"
        self.already_warned = set()

    def setup_parsers(self):
        # Each line contains data for a given station and period of time.
        com_parser = parser.define_parser(
            end_marker=lambda _l, _ln, _n: True,
            label=lambda line, _ln: line[0:4].isnumeric(),
            parser_def={
                True: {
                    "parser": self.parse_com_line,
                    "fields": {
                        "station": (0, 4),
                        "start_day": (5, 7),
                        "start_month": (8, 10),
                        "start_year": (11, 15),
                        "end_day": (16, 18),
                        "end_month": (19, 21),
                        "end_year": (22, 26),
                        "pulse_width": (27, 30),
                        "detector_type": (30, 35),
                        "op_mod": (35, 39),
                        "e_lev": (40, 44),
                        "dummy1": (44, 48),
                        "dummy2": (49, 52),
                        "lcmh": (53, 56),
                        "lcml": (57, 60),
                        "lcm": (61, 64),
                        "iflg": (65, 67),
                    },
                }
            },
        )
        return itertools.repeat(com_parser)

    def parse_com_line(self, line, _):
        """Parse one line from the center of mass-file

        Center of mass-data is stored in a dict with station CDPs as keys. We store the valid start and end-dates as
        proper Python datetimes.

        Args:
            line (Dict):   Data from one line in the center of mass-file.
        """
        line["start_time"] = datetime(
            int(line.pop("start_year")), int(line.pop("start_month")), int(line.pop("start_day")), 0, 0, 0
        )
        line["end_time"] = datetime(
            int(line.pop("end_year")), int(line.pop("end_month")), int(line.pop("end_day")), 23, 59, 59
        )
        self.data.setdefault(line.pop("station"), list()).append(line)

    def __missing__(self, key):
        """Look up the correct time interval for a given station

        We slightly abuse the __missing__-feature of dicts to implement lookup in date-intervals. If the lookup is
        given both a station CDP and a Time-object, we look up center of mass in the correct date interval.

        Args:
            key (Tuple):  Pair of station CDP (String) and time epoch (Time).

        Returns:
            Float:   Center of mass correction in meters.
        """
        try:
            station, time = key
        except ValueError:
            return self.data[key]

        com_info = self.data.get(station, [])
        for info_line in com_info:
            if info_line["start_time"] <= time.utc.datetime <= info_line["end_time"]:
                return float(info_line["lcm"]) * unit.mm2m

        if station not in self.already_warned:
            log.warn(
                "Missing center of mass data for CDP '{}'. Using default for {} ({} meters).",
                station,
                self.vars["satellite"].title(),
                DEFAULT_COM[self.vars["satellite"]],
            )
            self.already_warned.add(station)
        return DEFAULT_COM[self.vars["satellite"]]
