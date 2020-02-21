"""A parser for reading SLR prediction files

Description:
------------

Reads data from files in the CPF file format as defined in http://ilrs.gsfc.nasa.gov/docs/2006/cpf_1.01.pdf

"""

# Standard library imports
import itertools

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit
from midgard.parsers._parser_chain import ChainParser, ParserDef


@plugins.register
class SlrPredictionParser(ChainParser):
    """A parser for reading SLR prediction files (CPF format)
    """

    def setup_parser(self):
        # Data are organized in tables of positions (and sometimes velocities)
        header_parser = ParserDef(
            end_marker=lambda _l, _ln, nextline: nextline[0:2].isnumeric(),
            label=lambda line, _ln: line[0:2].upper(),
            parser_def={
                "H1": {
                    "parser": self.parse_default,
                    "fields": {
                        "cpf": (3, 6),
                        "format_version": (7, 9),
                        "ephemeris_source": (11, 14),
                        "year_of_ephemeris_production": (15, 19),
                        "month_of_ephemeris_production": (20, 22),
                    },
                },
                "H2": {"parser": self.parse_default, "fields": {"satellite_id": (3, 11)}},
            },
        )
        # Each line contains information about the satellite at a given time.
        orbit_parser = ParserDef(
            end_marker=lambda _l, _ln, nextline: nextline[0:2] == "99",
            label=lambda line, _ln: line[0:2],
            parser_def={
                "10": {
                    "parser": self.parse_position,
                    "fields": [
                        "record_type",
                        "direction_flag",
                        "mjd",
                        "seconds_of_day",
                        "leap_second_flag",
                        "pos_x",
                        "pos_y",
                        "pos_z",
                    ],
                },
                "20": {
                    "parser": self.parse_velocity,
                    "fields": ["record_type", "direction_flag", "vel_x", "vel_y", "vel_z"],
                },
            },
        )

        return itertools.chain([header_parser], itertools.repeat(orbit_parser))

    def parse_default(self, line, _):
        """Add the contents of line to data

        Args:
            line: Dict containing the fields of a line.
        """
        self.data.update(line)

    def parse_position(self, line, cache):
        mjd1 = int(line.pop("mjd"))
        mjd2 = float(line.pop("seconds_of_day")) * Unit.sec2day
        line["pos"] = np.array([float(line.pop("pos_x")), float(line.pop("pos_y")), float(line.pop("pos_z"))])

        self.data.setdefault("positions", dict()).setdefault((mjd1, mjd2), dict()).update(line)

    def parse_velocity(self, line, cache):
        """TODO: Add something here later. Some old files might contain velocities
        """
        pass
