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
    """A parser for reading a VLBI correlator report"""

    def __init__(self, file_path, encoding=None):
        super().__init__(file_path, encoding)
        self.raw = {}
        self.supported_version = -1

    def read_data(self):
        """Parse the *mk4.hist file

        self.data will be populated with information from the file if a supported version is found
        """
        self._find_version()
        with open(self.file_path, mode="rt") as fid:
            if self.supported_version == 2:
                self._parse_file_2(fid)
            elif self.supported_version == 3:
                self._parse_file_3(fid)

        if self.supported_version > 0:
           self._organize_data()

    def _find_version(self):
        """ Read the first 10 lines to guess format of content
        """
        with open(self.file_path, mode="r") as fid:
            line_count = 0
            for line in fid:
                line_count += 1
                if line.startswith("%CORRELATOR_REPORT_FORMAT"):
                    line = line.split()
                    self.supported_version = int(line[-1].strip())
                    # Found a specified version
                    return
                if line.startswith("+HEADER"):
                    # Found a version from before format is specified. Let's call this version 2
                    self.supported_version = 2
                    return
                if line_count > 10:
                    # Have not found anything useful. Give up
                    return

    def _parse_file_3(self, fid):
        """ Parse %CORRELATOR_REPORT_FORMAT 3
        """
        for line in fid:
            block = ""
            # Only read +'STATIONS and +QCODES. Skip everything else for now
            if not line.strip():
                # Skip empty lines
                continue
            if line.startswith("*"):
                # Skip comments
                continue
            if line.startswith("%"):
                # Skip format specifier
                continue
            if not line.startswith("+"):
                continue
                
            if line.startswith("+STATIONS") or line.startswith("+QCODES"):
                block = line[1:].strip()
                self.raw.setdefault(block, {})
                self._parse_block_3(fid, block)
            
                
    def _parse_block_3(self, fid, block):
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
        if empty.strip():
            # empty line was not empty. Assume this is header
            headers = empty.split()
        else:
            headers = fid.readline().split()
        for h in headers:
            self.raw[block].setdefault(h, [])
        dashes = fid.readline()
        for line in fid:
            if line.startswith("*"):
                continue
            # For some reports the dashed line contains spaces
            if not line.strip() or set(line.strip().replace(" ", "")) == {'-'}:
                return
            line = line.split()
            for h, l in zip(headers, line):
                try:
                    self.raw[block][h].append(int(l))
                except ValueError:
                    self.raw[block][h].append(l)

    def _parse_file_2(self, fid):
        """ Parse precursor to vlbi correlator report format.
            Search for +STATION_NOTES and +QCODES blocks
        """
        station_block = False
        for line in fid:
            if not line.strip():
                continue
            if line.startswith("*"):
                continue
            if line.startswith("+QCODES"):
                # Same format as QCODES block in version 3
                block = line[1:].strip()
                self.raw[block] = {}
                self._parse_block_3(fid, block)
                continue
            if line.startswith("+STATION_NOTES") or line.startswith("+STATION NOTES"):
                # Found station notes block
                station_block = True
                self.raw["STATIONS"] = {}
                self.raw["STATIONS"]["station"] = []
                self.raw["STATIONS"]["name"] = []
                self.raw["STATIONS"]["mk4"] = []
                continue
            if line.startswith("+"):
                # Start of any other block
                station_block = False
                continue
            if station_block:
                # Inside station block. Parse lines
                # Typical format: BADARY   (Bd/B): Ok.
                if line.find(":") > 0 and line.find("(") > 0 and line.find(")") > 0 and line.find("/") > 0 and not line.startswith("    "):
                    line, _ = line.split(":", maxsplit=1)
                    line = line.split("(")
                    self.raw["STATIONS"]["station"].append(line[1][0:2])
                    self.raw["STATIONS"]["name"].append(line[0].strip())
                    self.raw["STATIONS"]["mk4"].append(line[1][3:4])    


    def _organize_data(self):
        """ Copy content from self.raw to self.data and organize columns
        """
        if not self.raw["STATIONS"] or not self.raw["QCODES"]:
            # station or qcodes block where empty. Give up
            return
        self.data["stations"] = {}
        sta2letter = dict(zip(self.raw["STATIONS"]["station"], self.raw["STATIONS"]["mk4"]))
        name2letter = dict(zip(self.raw["STATIONS"]["name"], self.raw["STATIONS"]["mk4"]))
        self.data["stations"].update(sta2letter)
        self.data["stations"].update(name2letter)
        
        # Skip the 'total' column (Explains [:-1] indexing)
        keys = list(self.raw["QCODES"].keys())[:-1]
        keys0 = keys.pop(0)
        # For some correlator reports the keys[0] is different from "bl:band"
        bl = "bl"
        band = "band"
        data_type = [(k, "i4") for k in keys]
        data_type.insert(0, (band, "U2"))
        data_type.insert(0, (bl, "U2"))
        # Skip the 'total' row if present (will not be present if it there was a dashed line before the total row)
        l = len(self.raw["QCODES"][keys0])
        if self.raw["QCODES"][keys0][-1].lower().startswith("tot"):
            # Remove 'total' row
            l = -1
        else:
            l = len(self.raw["QCODES"][keys0])
        num_entries = len(self.raw["QCODES"][keys0][:l])
        self.data["qcodes"] = np.empty(num_entries, dtype=data_type)
        for k in keys:
            self.data["qcodes"][k] = self.raw["QCODES"][k][:l]
        split_bl_band = np.char.split(self.raw["QCODES"][keys0][:l], ":")
        bl_data, band_data = np.array(list(zip(*split_bl_band)))
        self.data["qcodes"][bl] = bl_data
        self.data["qcodes"][band] = band_data
