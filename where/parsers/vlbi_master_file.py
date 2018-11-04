"""A parser for reading the VLBI master file.

Description:
------------

Reads the VLBI master file which contains information about all the planned and executed sessions each year.

"""

# Standard library imports
import itertools

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import cache
from where.lib import config
from where.lib import log
from where.parsers._parser_chain import ParserDef, ChainParser


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
                    "delimiter": "|",
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
        self.data.setdefault((doy, session), dict()).update(line)
