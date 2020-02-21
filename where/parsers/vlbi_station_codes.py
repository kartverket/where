"""A parser for reading data from the VLBI station codes file

Example:
--------

    from where.parsers import VlbiStationCodesParser
    parser = VlbiStationCodesParser(file_key)
    parser.process_data()

Description:
------------

Reads data from the VLBI station codes file.

"""

# Standard library imports
import itertools

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_chain import ParserDef, ChainParser


@plugins.register
class VlbiStationCodesParser(ChainParser):
    """A parser for reading data from EOP files
    """

    def setup_parser(self):

        # Each line contains identifiers for a station
        station_parser = ParserDef(
            end_marker=lambda _l, _ln, _n: True,
            label=lambda line, _ln: not line.startswith("*"),
            parser_def={
                True: {
                    "parser": self.parse_station,
                    "fields": {
                        "ivscode": (0, 3),
                        "name": (4, 12),
                        "domes": (13, 18),
                        "marker": (18, 22),
                        "cdp": (23, 27),
                        "description": (28, 150),
                    },
                }
            },
        )

        return itertools.chain(itertools.repeat(station_parser))

    def parse_station(self, line, _):
        """Parse one line of station data

        Note: Some stations does not have a cdp number or a name.

        Args:
            line:  Dict containing the fields of a line.
        """
        line_copy = line.copy()
        name = line_copy.pop("name")
        cdp = line.pop("cdp")

        if cdp == "----":
            # Create a unique key
            cdp = "-" + line["ivscode"] + "-"
            line_copy["cdp"] = cdp

        if name == "--------":
            # Create a unique key
            name = "---" + line["ivscode"] + "---"
            line["name"] = name

        # If there are multiple entries with the same cdp number use the first
        if cdp not in self.data:
            # Store data twice with different keys
            self.data[cdp] = line
            self.data[name] = line_copy
