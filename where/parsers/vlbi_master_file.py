"""A parser for reading the VLBI master file.

Description:
------------

Reads the VLBI master file which contains information about all the planned and executed sessions each year.




"""

# Standard library imports
from datetime import date, timedelta
from collections import defaultdict
import itertools
import re

# Where imports
from where.lib import cache
from where.lib import config
from where.lib import log
from where.parsers._parser_chain import ParserDef, ChainParser
from where.lib import plugins


@plugins.register
class VlbiMasterFile(ChainParser):
    """A parser for reading VLBI antenna information
    """

    def setup_parser(self):
        # Each line contains identifiers for a session
        session_parser = ParserDef(
            end_marker=lambda _l, _ln, _n: True,
            label=lambda line, _ln: line.startswith("|"),
            parser_def={
                True: {
                    "parser": self.parse_session,
                    "fields": [
                        None,
                        "session_name",
                        "session_code",
                        None,
                        "doy",
                        "start",
                        "duration",
                        "stations",
                        "scheduler",
                        "correlator",
                        None,
                        None,
                        "session",
                        None,
                        None,
                        None,
                        None,
                    ],
                    "delimiter": "|",
                }
            },
        )

        return itertools.chain(itertools.repeat(session_parser))

    def parse_session(self, line, _):
        """Parse one line of antenna information

        Args:
            line:  Dict containing the fields of a line.
        """
        doy = int(line.pop("doy"))
        session = line.pop("session")
        reg_hits = re.search("\d", line["session_code"])
        num_idx = reg_hits.start() if reg_hits else len(line["session_code"])
        line["session_type"] = line["session_code"][:num_idx]
        self.data.setdefault((doy, session), dict()).update(line)

    @cache.function
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
        date, session = key
        log.warn("Session '{}' not found in master file for {}", session, date.strftime(config.FMT_date))

        return defaultdict(default_factory="")
