"""A parser for reading VLBI antenna information

Description:
------------

Reads data information about VLBI antenna dimensions, type and thermal expansion properties for the different antenna
construction elements.




"""

# Standard library imports
import itertools

# Where imports
from where.parsers._parser_chain import ParserDef, ChainParser
from where.lib import plugins


@plugins.register
class VlbiAntennaInfoParser(ChainParser):
    """A parser for reading VLBI antenna information
    """

    def setup_parser(self):
        # Each line contains identifiers for a station
        station_parser = ParserDef(
            end_marker=lambda _l, _ln, _n: True,
            label=lambda line, _ln: line.startswith("A"),
            parser_def={
                True: {
                    "parser": self.parse_station,
                    "fields": [
                        None,
                        "ivsname",
                        "focus",
                        "mount",
                        "radome",
                        "quality",
                        "reference_temperature",
                        "sin_amplitude",
                        "cos_amplitude",
                        "reference_pressure",
                        "diameter",
                        "height_foundation",
                        "depth_foundation",
                        "coefficient_foundation",
                        "fixed_axis",
                        "coefficient_fixed_axis",
                        "axis_offset",
                        "coefficient_axis_offset",
                        "distance_antenna_vertex",
                        "coefficient_distance",
                        "height_focus",
                        "coefficient_height_focus",
                    ],
                }
            },
        )

        return itertools.chain(itertools.repeat(station_parser))

    def parse_station(self, line, _):
        """Parse one line of antenna information

        Args:
            line:  Dict containing the fields of a line.
        """
        ivsname = line.pop("ivsname")
        self.data[ivsname] = dict()
        for key in line:
            try:
                self.data[ivsname][key] = float(line[key])
            except ValueError:
                self.data[ivsname][key] = line[key]
