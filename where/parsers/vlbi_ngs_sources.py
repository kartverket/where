"""A parser for reading VLBI source data from NGS files

Description:
------------

Reads data from files in the NGS file format as defined in http://lacerta.gsfc.nasa.gov/mk5/help/dbngs_format.txt
(revision date June 11, 2007) [1].

References:
-----------

[1] NGS file format.
    http://lacerta.gsfc.nasa.gov/mk5/help/dbngs_format.txt




$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""

# Standard library imports
import itertools

# External library imports
import numpy as np

# Where imports
from where.parsers._parser_chain import ParserDef, ChainParser
from where.lib import plugins
from where.lib.unit import unit


@plugins.register
class VlbiNgsSourcesParser(ChainParser):
    """A parser for reading VLBI Sites from the header of an NGS file
    """

    #
    # PARSERS for reading each line of the NGS file.
    #
    def setup_parser(self):
        # Ignore the header, first two lines
        header_parser = ParserDef(end_marker=lambda _l, line_num, _n: line_num == 2, label=None, parser_def=None)

        # Skip stations
        station_parser = ParserDef(end_marker=lambda line, _ln, _n: line == "$END", label=None, parser_def=None)

        # Each line defines a radio source
        source_parser = ParserDef(
            end_marker=lambda line, _ln, _n: line == "$END",
            label=lambda line, _ln: line != "$END",
            parser_def={
                True: {
                    "parser": self.parse_radio_source,
                    "fields": {
                        "name": (0, 8),
                        "ra_hrs": (10, 12),
                        "ra_mins": (13, 15),
                        "ra_secs": (16, 28),
                        "dec_degs": (29, 32),
                        "dec_mins": (33, 35),
                        "dec_secs": (36, 48),
                    },
                }
            },
        )

        return itertools.chain([header_parser, station_parser, source_parser])

    def parse_radio_source(self, line, _):
        """Read station position

        Reads the station position from the NGS file.

        Args:
            line:  Input data from NGS file
        """
        src_name = line["name"].replace(".", "dot")
        self.data[src_name] = dict()
        self.data[src_name]["ra"] = unit.hms_to_rad(
            float(line["ra_hrs"]), int(line["ra_mins"]), float(line["ra_secs"])
        )
        self.data[src_name]["dec"] = unit.dms_to_rad(
            float(line["dec_degs"].replace(" ", "")), int(line["dec_mins"]), float(line["dec_secs"])
        )
