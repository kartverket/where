"""A parser for reading VLBI correlation report within a vgosDb

Description:
------------

Reads data from files in the vgosDb files as defined in [1]. The data is organized in multiple smaller database
files based on netCDF.

References:
-----------

..[1] vgosDb format
    ftp://gemini.gsfc.nasa.gov/pub/misc/jmg/VLBI_Structure_2013Jun11.pdf

    new url needed

"""

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser import Parser
from midgard import parsers

# Where imports
from where.lib import log


@plugins.register
class CorrelationReportParser(Parser):
    """A parser for reading VLBI data from a VGOS database"""

    def __init__(self, file_path, encoding=None):
        super().__init__(file_path, encoding)
        self.raw = {}

    def read_data(self):
        """Parse the vgosdb wrapper file

        self.data will be populated with information from the netcdf files
        """
        with open(self.file_path, mode="rt") as fid:
            self._parse_file(fid)

        self._organize_data()

    def _parse_file(self, fid):
        for line in fid:
            block = ""
            # Only read +'STATIONS and +QCODES. Skip everything else for now
            if not line.strip():
                # Skip empty lines
                continue
            if line.startswith("*"):
                # Skip comments
                continue
            if line.startswith("%CORRELATION_REPORT_FORMAT"):
                line = line.split()
                if line[-1] != 3:
                    log.warn(f"Unknown correlation report format. Version {line[-1]}")
            if not line.startswith("+"):
                continue
                
            if line.startswith("+STATIONS") or line.startswith("+QCODES"):
                block = line[1:].strip()
                self.raw.setdefault(block, {})
                self._parse_block(fid, block)
            
                
    def _parse_block(self, fid, block):
        """ Parse a data block.
        
        Expected format:
        +<block>
        <empty line>
        <column headers>
        -----------------
        <columns with data>
        <empty line> # stop parsing at this point
        <comments>
        <empty line>
        """
        empty = fid.readline()
        headers = fid.readline().split()
        for h in headers:
            self.raw[block].setdefault(h, [])
        dashes = fid.readline()
        for line in fid:
            if not line.strip():
                return
            line = line.split()
            for h, l in zip(headers, line):
                try:
                    self.raw[block][h].append(int(l))
                except ValueError:
                    self.raw[block][h].append(l)


    def _organize_data(self):
        """ Copy content from self.raw to self.data and organize columns
        """
        self.data["stations"] = {}
        sta2letter = dict(zip(self.raw["STATIONS"]["station"], self.raw["STATIONS"]["mk4"]))
        name2letter = dict(zip(self.raw["STATIONS"]["name"], self.raw["STATIONS"]["mk4"]))
        self.data["stations"].update(sta2letter)
        self.data["stations"].update(name2letter)
        
        self.data["qcodes"] = {}
        # Skip the 'total' column and 'total' row
        keys = list(self.raw["QCODES"].keys())[:-1]
        data_type = [(k, "i4") for k in keys[1:]]
        data_type.insert(0, (keys[0], "U5"))
        num_entries = len(self.raw["QCODES"][keys[0]][:-1])
        self.data["qcodes"] = np.empty(num_entries, dtype=data_type)
        for k in keys:
            self.data["qcodes"][k] = self.raw["QCODES"][k][:-1]

