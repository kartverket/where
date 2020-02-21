"""A parser for reading VLBI eccentricity vectors from file

Description:
------------

Reads the VLBI eccentricity vector from file.

"""

# Standard library imports
from datetime import datetime
import re

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_line import LineParser


@plugins.register
class VlbiEccentricityParser(LineParser):
    """A parser for reading VLBI eccentricity vectors from file
    """

    def setup_parser(self):
        def str2date(s):
            return datetime.strptime(re.sub("[^0-9]", " ", s.decode()), "%Y %m %d %H %M")

        def empty(s):
            return s.decode().replace("?", "")

        return dict(
            names="name, site_id, start, end, v1, v2, v3, coord_type",
            dtype=("U10", "U5", object, object, "f8", "f8", "f8", "U6"),
            delimiter=(10, 5, 18, 18, 12, 12, 12, 6),
            skip_header=1,
            skip_footer=1,
            comments="$",
            converters={0: empty, 2: str2date, 3: str2date},
            autostrip=True,
        )

    def structure_data(self):
        for item in self._array:
            key = (str(item["name"]).strip().replace(" ", "_"), str(item["site_id"]).strip())

            if item["coord_type"] == "NEU":
                # Swap NEU to ENU
                self.data.setdefault(key, {}).setdefault((item["start"], item["end"]), {}).update(
                    dict(vector=(item["v2"], item["v1"], item["v3"]), coord_type="ENU")
                )
            else:
                self.data.setdefault(key, {}).setdefault((item["start"], item["end"]), {}).update(
                    dict(vector=(item["v1"], item["v2"], item["v3"]), coord_type=item["coord_type"])
                )
            self.data[key]["name"] = item["name"]
