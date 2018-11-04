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

# Where imports
from where.parsers._parser_chain import ParserDef, ChainParser
from where.lib import cache
from where.lib import log


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

        Args:
            line:  Dict containing the fields of a line.
        """
        line_copy = line.copy()
        name = line_copy.pop("name")
        cdp = line.pop("cdp")
        # If there are multiple entries with the same cdp number use the first
        if cdp not in self.data:
            # Store data twice with different keys
            self.data[cdp] = line
            self.data[name] = line_copy

    @cache.function
    def __missing__(self, cdp):
        """Handle missing keys

        Give a warning and return a consistent value when eccentricity vectors are missing from the file.

        The special __missing__ method is called when doing a __getitem__ (i.e. `dict_name[key]`) lookup on a
        dictionary and the key is missing.

        The caching is used mainly as a simple way of only warning about missing eccentricity data. This could also
        have been implemented by keeping a class set containing which CDPs we have already warned about.

        Args:
            cdp (String):  cdp for the missing station

        Returns:
            Dict:  Dummy information
        """
        log.warn("Missing station codes for '{}'.", cdp)
        return dict(ivscode="", ivsname="", domes="", marker="", description="")
