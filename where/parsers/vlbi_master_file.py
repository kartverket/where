"""A parser for reading the VLBI master file.

Description:
------------

Reads the VLBI master file which contains information about all the planned and executed sessions each year.

"""

# Standard library imports
from datetime import datetime
import itertools

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_chain import ParserDef, ChainParser


@plugins.register
class VlbiMasterFile(ChainParser):
    """A parser for reading VLBI antenna information
    """

    def __init__(self, file_path, encoding=None):
        # Read first line of file to determine the version of the file format
        with open(file_path, mode="rt", encoding=encoding) as fid:
            # Typical header line:
            # ## Master file format version 1.0           2001.08.21 CCT&NRV
            header = fid.read().split()
            self.version = header[5]

        super().__init__(file_path, encoding)

    def setup_parser(self):

        if self.version == "1.0":
            # Each line contrains identifiers for a session
            session_parser = ParserDef(
                end_marker=lambda _l, _ln, _n: True,
                label=lambda line, _ln: line.startswith("|"),
                parser_def={
                    True: {
                        "parser": self.parse_session,
                        "delimiter": "[|]", # Regex syntax
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
                            "status",
                            None,
                            "dbc",
                            None,
                            None,
                            None,
                            None,
                        ],
                    }
                },
            )
        elif self.version == "2.0":
            # TODO: new version of master file
            # Each line contrains identifiers for a session
            session_parser = ParserDef(
                end_marker=lambda _l, _ln, _n: True,
                label=lambda line, _ln: line.startswith("|"),
                parser_def={
                    True: {
                        "parser": self.parse_session,
                        "delimiter": "[|]",
                        "fields": [
                            None,
                            "session_type",
                            "date",
                            "session_code",
                            "doy",
                            "start",
                            "duration",
                            "stations",
                            "scheduler",
                            "correlator",
                            "status",
                            "dbc",
                            "submitter",
                            "delay",
                        ],
                    }
                },
            )
        else:
            print(f"Unknown version of master file format: {self.version}")
            # Error finding version.

        return itertools.chain(itertools.repeat(session_parser))

    def parse_session(self, line, _):
        """Parse one line of antenna information

        Args:
            line:  Dict containing the fields of a line.
        """
        #doy = int(line.pop("doy"))
        session_code = line.pop("session_code")
        self.data.setdefault(session_code, dict()).update(line)
        self.data[session_code]["doy"] = int(self.data[session_code]["doy"]) 

        # Convert status information to a date object if possible
        status = line["status"].strip()
        if self.version == "1.0":
            try:
                # Avoid using %b, %-m and %-d in strptime. Platform specific and weird locale behaviour
                months = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
                yy = int(status[0:2])
                century = 2000 if yy < 50 else 1900
                yyyy = str(yy + century)
                MMM = status[2:5]
                dd = status[5:7]
                mm = str(months.index(MMM) + 1).zfill(2)
                status_date = datetime.strptime(f"{yyyy}{mm}{dd}", "%Y%m%d").date()
                self.data[session_code]["status"] = status_date
            except (IndexError, ValueError):
                pass

        if self.version == "2.0":
            try:
                status_date = datetime.strptime(status, "%Y%m%d").date()
                self.data[session_code]["status"] = status_date
            except ValueError:
                pass
