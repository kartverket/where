"""Get current master schedule for VLBI

Description:
------------


"""
from collections import UserDict, defaultdict

# Where imports
from where import parsers
from where.lib import config
from where.lib import log
from where.lib import plugins
from where.lib import cache


@plugins.register
def get_vlbi_master_schedule(rundate=None):
    """Read master file

    If rundate is not specified, used file_vars that are already set.
    """

    file_vars = None if rundate is None else config.date_vars(rundate)
    parser = parsers.parse_key("vlbi_master_file", file_vars=file_vars)

    return VlbiMasterSchedule(parser.as_dict(), file_vars["yyyy"])


class VlbiMasterSchedule(UserDict):
    def __init__(self, data, year):
        self.data = data
        self.year = year

    def __missing__(self, key):
        """Handle missing keys

        Give a warning and return a consistent value when session is missing in master file.

        The special __missing__ method is called when doing a __getitem__ (i.e. `dict_name[key]`) lookup on a
        dictionary and the key is missing.

        Args:
            key (String):  key for the missing value.

        Returns:
            Dict:  Dummy information, blank strings for all fields.
        """
        doy, session = key
        log.warn("Session '{}' not found in master file for {} doy {}", session, self.year, doy)

        return defaultdict(default_factory="")

    def list_sessions(self, date):
        """List sessions available at a given date

        Args:
            date (date):  The given date.

        Returns:
            List of strings:  Names of sessions available at a given date.
        """
        doy = date.timetuple().tm_yday
        return [s for d, s in self.data.keys() if d == doy]
