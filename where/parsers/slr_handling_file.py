"""A parser for reading SLR handling file

Description:
------------

Asdf

References:
-----------

http://ilrs.dgfi.tum.de/fileadmin/data_handling/ILRS_Data_Handling_File.snx



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""
# Standard library imports
from datetime import datetime, timedelta
import itertools

# External library imports

# Where imports
from where.parsers import parser
from where.lib import plugins


@plugins.register
class SlrHandlingFileParser(parser.ParserDict):
    """A parser for reading SLR handling file
    """

    def __init__(self, time):
        super().__init__()

        self.start_time = min(time.utc).datetime
        self.end_time = max(time.utc).datetime

    def setup_parsers(self):
        # Each line contains data for a given station and period of time.
        handling_parser = parser.define_parser(
            end_marker=lambda _l, _ln, _n: True,
            label=lambda line, _ln: line[1:5].isnumeric() and line[18:19].isnumeric(),
            parser_def={
                True: {
                    "parser": self.parse_handling_line,
                    "fields": {
                        "station": (1, 5),
                        "unit": (10, 12),
                        "start_time": (17, 29),
                        "end_time": (30, 42),
                        "code": (43, 44),
                    },
                }
            },
        )
        return itertools.repeat(handling_parser)

    def parse_handling_line(self, line, _):
        # The infinity mark on file is 00:000:00000
        if line["start_time"] == "00:000:00000":
            start_time = datetime.min
        else:
            start_time = (
                datetime.strptime(line["start_time"][:6], "%y:%j") + timedelta(seconds=int(line["start_time"][7:]))
            )
        if line["end_time"] == "00:000:00000":
            end_time = datetime.max
        else:
            end_time = (
                datetime.strptime(line["end_time"][:6], "%y:%j") + timedelta(seconds=int(line["end_time"][7:]))
            )

        interval = (start_time, end_time)
        if (start_time <= self.start_time <= end_time) or (self.start_time <= start_time <= self.end_time):
            self.data.setdefault(line["station"], {}).setdefault(line["code"], []).append((interval, line["unit"]))
