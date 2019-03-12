"""
Description:

    Reads station-dependent information from file
    which contains informaton on which data sources
    are reliable or not.

References:

    http://ilrs.dgfi.tum.de/fileadmin/data_handling/ILRS_Data_Handling_File.snx

"""

# Standard library imports
from collections import UserDict

# Where imports
from where import parsers
from where.lib import plugins


@plugins.register
def get_handling_file(time):
    """Read station-dependent info
    """
    handling_parser = parsers.parse_key("slr_handling_file")
    return HandlingFile(handling_parser.as_dict(), time)


class HandlingFile(UserDict):
    def __init__(self, data, time):
        super().__init__()
        self.start_time = min(time.utc).datetime
        self.end_time = max(time.utc).datetime
        self._pick_data(data)

    def _pick_data(self, data):
        for site, site_data in data.items():
            for code, intervals in site_data.items():
                for interval, info in intervals:
                    start_time, end_time = interval
                    if (start_time <= self.start_time <= end_time) or (self.start_time <= start_time <= self.end_time):
                        self.data.setdefault(site, {}).setdefault(code, []).append((interval, info))
