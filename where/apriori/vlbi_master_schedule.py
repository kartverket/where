"""Get current master schedule for VLBI

Description:
------------


"""
import re
from collections import UserDict, defaultdict

# Where imports
from where import parsers
from where.lib import config
from where.lib import log
from where.lib import plugins


@plugins.register
def get_vlbi_master_schedule(rundate=None):
    """Read master file

    If rundate is not specified, used file_vars that are already set.
    """
    if rundate:
        file_vars = config.date_vars(rundate)
        parser = parsers.parse_key("vlbi_master_file", file_vars=file_vars)
        return VlbiMasterSchedule(parser.as_dict(), file_vars["yyyy"])

    else:
        return VlbiMasterSchedule


class VlbiMasterSchedule(UserDict):
    def __init__(self, data, year):
        super().__init__(data)
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
        log.warn(f"Session {session!r} not found in master file for {self.year} doy {doy}")

        return defaultdict(default_factory="")

    def list_sessions(self, date, session_types=None):
        """List sessions available at a given date

        Args:
            date (date):           The given date.
            session_types (List):  Optional filter, only show sessions in this list.

        Returns:
            List of strings:  Names of sessions available at a given date.
        """
        doy = date.timetuple().tm_yday
        if session_types:
            # Filter on session types in addition to date
            return [
                s
                for ((d, s), v) in self.data.items()
                if d == doy and self.session_type(v["session_code"]) in session_types
            ]
        return [s for d, s in self.data.keys() if d == doy]

    @staticmethod
    def session_type(session_code):
        reg_hits = re.search("\d", session_code)
        num_idx = reg_hits.start() if reg_hits else len(session_code)
        return session_code[:num_idx]
