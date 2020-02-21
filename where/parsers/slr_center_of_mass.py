"""A parser for reading SLR center of mass corrections from file

Description:
------------

Asdf

"""

# Standard library imports
from datetime import datetime, time

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_line import LineParser


@plugins.register
class SlrComParser(LineParser):
    """A parser for reading SLR center of mass corrections from files
    """

    def setup_parser(self):
        def str2date_min(s):
            return datetime.strptime(s.decode(), "%d %m %Y")

        def str2date_max(s):
            return datetime.combine(datetime.strptime(s.decode(), "%d %m %Y"), time.max)

        return dict(
            names=[
                "site_code",
                "start",
                "end",
                "pulse_width",
                "detector_type",
                "op_mod",
                "e_lev",
                "dummy1",
                "dummy2",
                "lcmh",
                "lcml",
                "lcm",
                "iflg",
            ],
            delimiter=(4, 11, 11, 4, 5, 5, 4, 4, 4, 4, 4, 4, 3),
            dtype=("U4", object, object, "i8", "U5", "U5", "f8", "i8", "i8", "i8", "i8", "i8", "i8"),
            converters={1: str2date_min, 2: str2date_max},
            autostrip=True,
        )

    def structure_data(self):
        fields = self._array.dtype.names[1:]
        for item in self._array:
            self.data.setdefault(item["site_code"], list()).append({f: item[f] for f in fields})
