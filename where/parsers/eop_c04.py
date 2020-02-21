"""A parser for reading data from EOP files

Description:
------------

Reads data from EOP files.

"""

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_line import LineParser


@plugins.register
class EopC04Parser(LineParser):
    """A parser for reading data from EOP files
    """

    def setup_parser(self):
        # TODO: better solution
        skip_header = 14 if self.file_path.name.endswith("now") else 7
        return dict(
            names=["year", "month", "day", "mjd", "x", "y", "ut1_utc", "lod", "dx", "dy"],
            skip_header=skip_header,
            usecols=(3, 4, 5, 6, 7, 8, 9),
        )

    def structure_data(self):
        self.data = {
            item["mjd"]: dict(
                x=item["x"],
                y=item["y"],
                ut1_utc=item["ut1_utc"],
                lod=item["lod"],
                dx=item["dx"],
                dy=item["dy"],
                source="c04",
            )
            for item in self._array
        }
