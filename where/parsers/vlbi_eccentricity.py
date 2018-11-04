"""A parser for reading VLBI eccentricity vectors from file

Description:
------------

Reads the VLBI eccentricity vector from file.

"""

# Standard library imports
from datetime import datetime, time
import itertools
import re

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import cache
from where.lib import log
from where.parsers import parser


@plugins.register
class VlbiEccentricityParser(parser.ParserDict):
    """A parser for reading VLBI eccentricity vectors from file
    """

    #
    # PARSER for reading each line of the EOP file.
    #
    def setup_parsers(self):
        # Each line contains identifiers for a station
        station_parser = parser.define_parser(
            end_marker=lambda _l, _ln, _n: True,
            label=lambda line, ln: not (line.startswith("$") or line.startswith("#")),
            parser_def={
                True: {
                    "parser": self.parse_station,
                    "fields": {
                        "name": (2, 10),
                        "site_id": (11, 16),
                        "start": (17, 34),
                        "end": (35, 52),
                        "v1": (53, 64),
                        "v2": (64, 75),
                        "v3": (75, 85),
                        "type": (87, 91),
                    },
                }
            },
        )

        return itertools.chain(itertools.repeat(station_parser))

    def parse_station(self, line, _):
        """Parse one line of eccentricity information

        Only stores the eccentricities that are relevant according to self.rundate.

        Args:
            line (Dict):  The fields of a line.
        """
        site_id = line.pop("site_id")

        # Replace all non-numeric characters in dates with spaces (since the format is inconsistent)
        start = datetime.strptime(re.sub("[^0-9]", " ", line.pop("start")), "%Y %m %d %H %M")
        end = datetime.strptime(re.sub("[^0-9]", " ", line.pop("end")), "%Y %m %d %H %M")
        if datetime.combine(self.rundate, time.max) > start and datetime.combine(self.rundate, time.min) < end:
            self.data[site_id] = dict(
                start=start,
                end=end,
                vector=np.array((float(line.pop("v1")), float(line.pop("v2")), float(line.pop("v3")))),
            )
            self.data[site_id].update(line)

    @cache.function
    def __missing__(self, site_id):
        """Handle missing keys

        Give a warning and return a consistent value when eccentricity vectors are missing from the file.

        The special __missing__ method is called when doing a __getitem__ (i.e. `ecc[key]`) lookup on a dictionary and
        the key is missing.

        The caching is used mainly as a simple way of only warning about missing eccentricity data. This could also
        have been implemented by keeping a class set containing which CDPs we have already warned about.

        Args:
            site_id (String):  CDP-number for the missing eccentricity vector

        Returns:
            Dict:  Dummy information including a zero eccentricity vector
        """
        log.warn("Missing eccentricity data for site id '{}'. Vector set to zero.", site_id)
        return dict(vector=np.zeros(3), name="Name missing", type="XYZ", start=datetime.min, end=datetime.max)
