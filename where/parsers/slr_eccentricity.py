"""A parser for reading SLR eccentricity vectors from file

Description:
------------

Reads the SLR eccentricity vector from file.



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

# Standard library imports
import itertools
from datetime import datetime, timedelta, time

# External library imports
import numpy as np

# Where imports
from where.lib import cache
from where.lib import log
from where.parsers import parser
from where.lib import plugins


@plugins.register
class SlrEccentricityParser(parser.ParserDict):
    """A parser for reading SLR eccentricity vectors from file
    """

    def __init__(self, rundate):
        super().__init__(rundate)
        self.file_key = "eccentricity"

    #
    # PARSER for reading each line of the eccentricity file.
    #
    def setup_parsers(self):
        # Each line contains identifiers for a station
        station_parser = parser.define_parser(
            end_marker=lambda _l, _ln, _n: True,
            label=lambda line, ln: (
                not (line.startswith("$") or line.startswith("#") or line.startswith("*")) and line[42:45] == "XYZ"
            ),
            parser_def={
                True: {
                    "parser": self.parse_station,
                    "fields": {
                        "site_id": (1, 5),
                        "start": (16, 28),
                        "end": (29, 41),
                        "type": (42, 45),
                        "v1": (45, 54),
                        "v2": (54, 63),
                        "v3": (63, 72),
                    },
                }
            },
        )

        return itertools.chain(itertools.repeat(station_parser))

    def parse_station(self, line, _):
        """Parse one line of eccentricity information

        Only saves the eccentricity relevant according to self.rundate

        Args:
            line (Dict):  The fields of a line.
        """
        site_id = line.pop("site_id")
        start = line.pop("start").replace(":", " ")
        if start[3:6] == "000":
            start = "70 001 00000"

        end = line.pop("end").replace(":", " ")
        if end[3:6] == "000":
            end = "50 001 00000"

        start = datetime.strptime(start[:6], "%y %j") + timedelta(seconds=int(start[7:]))
        end = datetime.strptime(end[:6], "%y %j") + timedelta(seconds=int(end[7:]))

        if datetime.combine(self.rundate, time.max) > start and datetime.combine(self.rundate, time.min) < end:
            self.data[site_id] = dict(
                vector=np.array((float(line.pop("v1")), float(line.pop("v2")), float(line.pop("v3")))),
                start=start,
                end=end,
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
