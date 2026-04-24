"""Get current master schedule for VLBI

Description:
------------


"""
import re
from collections import UserDict, defaultdict
from datetime import date

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
        session_code = key
        log.warn(f"Session {session_code!r} not found in master file for {self.year}")

        return defaultdict(default_factory="")

    def list_sessions(self, date, session_types=None, skip_intensives=True):
        """List sessions available at a given date

        Args:
            date (date):           The given date.
            session_types (List):  Optional filter, only show sessions in this list.

        Returns:
            List of strings:  Names of sessions available at a given date.
        """
        doy = date.timetuple().tm_yday
        sessions = list()
        for sc, v in self.data.items():
            if v["doy"] == doy:
                # Filter on session types in addition to date
                # Convert all given session types to upper case to make a case insensitive comparison
                if session_types:
                    session_types = [st.upper() for st in session_types]
                    if self.where_session_type(sc).upper() not in session_types:
                        continue
                # Before 1992 intensive sessions and 24 hour sessions are listed in the same master file
                if skip_intensives and v["session_type"] == "INTENSIVE":
                    continue
                sessions.append(sc)
        return sessions

    def ready(self, session_code):
        """Returns True if a session is marked as submitted in the master file

        Also returns True if information is missing. A session may be ready be processed even if
        it is not ready. Example: A preliminary database is released while waiting for data from 
        a station.
        """
        status = self.data[session_code]["status"]
        is_date = isinstance(status, date)
        if is_date:
            # If the status field is a date the session is ready
            return is_date

        is_str = isinstance(status, str)
        if is_str:
            # Return True for empty string (information is missing), otherwise False
            return not bool(status.strip())
        return False

    def status(self, session_code):
        return self.data[session_code]["status"]

    @staticmethod
    def where_session_type(session_code):
        """Returns custom defintion of session type for a session.

        The custom definition of session type is different from the official definition of session type
        available in version 2 of the master file format. The custom definition is needed as long as offical
        historical master files are not converted to version 2.

        The custom session type is defined as as all letter in the session code until the first digit is
        encountered. The letters are converted to uppercase.

        Example: Session code: R11000, Session type: R

        Args:
            session_code     Session code as defined in the master file

        Returns:
            The custom defined session type
        """
        session_code = session_code.upper()
        reg_hits = re.search(r"\d", session_code)
        num_idx = reg_hits.start() if reg_hits else len(session_code)
        return session_code[:num_idx]
