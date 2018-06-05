"""A parser for reading data from EOP files

Description:
------------

Reads data from EOP files.



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

# Standard library imports
import itertools

# Where imports
from where.parsers._parser_chain import ParserDef, ChainParser
from where.lib import plugins


@plugins.register
class EopParser(ChainParser):
    """A parser for reading data from EOP files
    """

    def setup_parser(self):
        # Each line contains data for a given date
        date_parser = ParserDef(
            end_marker=lambda _l, _ln, _n: True,
            label=lambda line, _ln: line[0:4].isnumeric(),
            parser_def={
                True: {
                    "parser": self.parse_eop_line,
                    "fields": [
                        None, None, None, "mjd", "x", "y", "ut1_utc", "lod", "dx", "dy"  # ignore date, use MJD
                    ],
                }
            },
        )
        return itertools.repeat(date_parser)

    def parse_eop_line(self, line, _):
        """Parse one line of EOP data

        First identify the date, and then convert the rest of the fields to float values.

        Args:
            line:  Dict containing the fields of a line.
        """
        obsdate = int(line.pop("mjd"))
        self.data[obsdate] = {k: float(v) for k, v in line.items()}
