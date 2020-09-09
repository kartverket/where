"""A parser for reading VLBI axis offsets

Description:
------------

Reads data information about VLBI axis offsets provided by GSFC.

This file does not contain information about mount type

"""

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_line import LineParser


@plugins.register
class VlbiAxisOffsetParser(LineParser):
    """A parser for reading VLBI axis offset"""

    def setup_parser(self):
        return dict(
            names="ivsname, apriori_offset, aposteriori_offset, plussminus, ,ferr, comment",
            dtype=("U8", "f8", "f8", "U3", "f8", "U40"),
            delimiter=(8, 10, 10, 3, 8, 40),
            skip_header=1,
            skip_footer=1,
            comments="#",
            autostrip=True,
        )

    def structure_data(self):
        for item in self._array:
            key = str(item["ivsname"]).strip().replace(" ", "_")
            value = item["aposteriori_offset"]
            self.data[key] = value
