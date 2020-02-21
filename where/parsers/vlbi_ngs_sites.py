"""A parser for reading VLBI station data from NGS files

Description:
------------

Reads data from files in the NGS file format as defined in http://lacerta.gsfc.nasa.gov/mk5/help/dbngs_format.txt
(revision date June 11, 2007) [1].

References:
-----------

[1] NGS file format.
    http://lacerta.gsfc.nasa.gov/mk5/help/dbngs_format.txt

"""

# Standard library imports
import itertools

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_chain import ParserDef, ChainParser


@plugins.register
class VlbiNgsSitesParser(ChainParser):
    """A parser for reading VLBI Sites from the header of an NGS file
    """

    #
    # PARSERS for reading each line of the NGS file.
    #
    def setup_parser(self):
        # Ignore the header, first two lines
        header_parser = ParserDef(end_marker=lambda _l, line_num, _n: line_num == 2, label=None, parser_def=None)

        # Each line defines a station (site)
        station_parser = ParserDef(
            end_marker=lambda line, _ln, _n: line == "$END",
            label=lambda line, _ln: line != "$END" and "station",
            parser_def={
                "station": {
                    "parser": self.parse_station,
                    "fields": {"name": (0, 8), "pos_x": (10, 25), "pos_y": (25, 40), "pos_z": (40, 55)},
                }
            },
        )

        return itertools.chain([header_parser, station_parser])

    def parse_station(self, line, _):
        """Read station position

        Reads the station position from the NGS file.

        Args:
            line:  Input data from NGS file
        """
        self.data[line["name"]] = np.array([float(line["pos_x"]), float(line["pos_y"]), float(line["pos_z"])])
