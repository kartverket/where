"""A parser for VLBI Schedule files

Description:
------------

Reads a skd file to extract the number of scheduled observations per station. 

All other information in the skd file is ignored.

Does not read VEX format.

"""

# Standard library imports
import itertools

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser import Parser

# Some stations do not use their standard names in the sked files.
STATION_NAME_MAP = {
    "TIGO": "TIGOCONC",
    "CTVA": "CTVASTJ",
    }


@plugins.register
class VlbiSkdParser(Parser):
    """A parser for reading VLBI schedule files produced by SKED."""

    def read_data(self):
        """Reads the schedule file and stores the information in self.data"""
        with open(self.file_path, mode="rt") as fid:
            all_lines = self._parse_file(fid)
        
        sta_codes = self._parse_stations(all_lines)
        num_obs = self._parse_scans(all_lines)
        self.data = {STATION_NAME_MAP.get(sta_codes[k], sta_codes[k]):v for k, v in num_obs.items() if k in sta_codes}
        if not self.data:
            self.data_available = False

    def _parse_file(self, fid):
        """Store each line in the file in a dictionary based on blocks

        Each block is defined by with the "$" character and the name of the block.

        Returns:
            dict:        key: $<block_name>, value: list of lines in block
        """
        all_lines = dict()
        first_line = fid.readline()
        if not first_line.startswith("$EXPER"):
            # Content is probably VEX-format, which is not implemented
            return all_lines

        for line in fid:
            if not line:
                continue

            if line.startswith("*"):
                # Assume the line is a comment
                continue

            if line.startswith("$"):
                key = line.strip()
                all_lines.setdefault(key, [])
                continue

            all_lines[key].append(line)
        return all_lines    

    def _parse_stations(self, all_lines):
        """Parses the station block

        Args:
            all_lines (dict):    all blocks from file
        Returns:
            dict:    key: one letter station code, value: station name
        """ 
        sta_codes = {}
        # Station block may have different names
        lines = all_lines.get("$STATION") or all_lines.get("$STATIONS") or []
        for line in lines:
            if not line:
                continue
            if line.startswith("A "):
                line = line.split()
                sta_codes[line[1]] = line[2]
        return sta_codes

    def _parse_scans(self, all_lines):
        """Parses the scans block

        Args:
            all_lines (dict):    all blocks from file
        Returns:
            dict:    key: one letter station code, value: number of scheduled observations
        """
        num_obs = {}
        for line in all_lines.get("$SKED", []):
            line = line.split()
            if not line:
                continue
            # Every second letter indicates refers to a station that is in the scan
            stations = line[9][::2]
            for s in stations:
                num_obs.setdefault(s, 0)
                num_obs[s] += len(stations) - 1

        return num_obs