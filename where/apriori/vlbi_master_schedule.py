"""Get current master schedule for VLBI

Description:
------------


"""
import re
from collections import UserDict, defaultdict
from datetime import datetime

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import parsers
from where.lib import config
from where.lib import log


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

    def ready(self, date, session):
        # Avoid using %b, %-m and %-d in strptime. Platform specific and weird locale behaviour
        months = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
        key = (date.timetuple().tm_yday, session)
        status = self.data[key]["status"].strip()
        if not status:
            # Status information might be missing
            return True
        try:
            yy = int(status[0:2])
            century = 2000 if yy < 50 else 1900
            yyyy = str(yy + century)
            MMM = status[2:5]
            dd = status[5:7]
            mm = str(months.index(MMM) + 1).zfill(2)
            datetime.strptime(f"{yyyy}{mm}{dd}", "%Y%m%d")
            # Assumption: successful parsing of status field as date means it's ready for processing
            return True
        except (IndexError, ValueError):
            return False

    def status(self, date, session):
        key = (date.timetuple().tm_yday, session)
        return self.data[key]["status"]

    @staticmethod
    def session_type(session_code):
        reg_hits = re.search("\d", session_code)
        num_idx = reg_hits.start() if reg_hits else len(session_code)
        return session_code[:num_idx]
