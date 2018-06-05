"""A parser for reading VLBI data from vgosDb files

Description:
------------

Reads data from files in the vgosDb files as defined in [1]. The data is organized in multiple smaller database
files based on netCDF.

References:
-----------

..[1] vgosDb format
    ftp://gemini.gsfc.nasa.gov/pub/misc/jmg/VLBI_Structure_2013Jun11.pdf


Authors:
--------

* Ann-Silje Kirkvik <ann-silje.kirkvik@kartverket.no>

$Revision: 15235 $
$Date: 2018-06-01 09:45:02 +0200 (Fri, 01 Jun 2018) $
$LastChangedBy: kirann $
"""
# Standard library imports
import os

# External library imports
import numpy as np
from scipy import interpolate

# Where imports
from where import apriori
from where.lib import config
from where.lib import constant
from where.lib import files
from where.lib import log
from where.lib.time import Time
from where.parsers._parser import Parser
from where import parsers
from where.lib import plugins
from where.ext import sofa


@plugins.register
class VgosDbParser(Parser):
    """A parser for reading ocean tidal loading coefficients from BLQ-files
    """

    SKIP_BLOCKS_DEFAULT = ["History", "Process"]

    def read_data(self):
        if True:
            self.SKIP_BLOCKS = self.SKIP_BLOCKS_DEFAULT
        else:
            self.SKIP_BLOCKS = []

        with files.open_path(self.file_path, mode="rt") as fid:
            self._parse_file(fid)

    def _parse_file(self, fid):
        for line in fid:
            if not line or line.startswith("!"):
                continue
            line = line.split()
            if "Begin" in line[0] and line[1] not in self.SKIP_BLOCKS:
                self._parse_block(fid, line[1], name=" ".join(line[2:]))

    def _parse_block(self, fid, block, name=""):
        print("Parsing {} {}".format(block, name))
        directory = ""
        for line in fid:
            if not line or line.startswith("!"):
                continue
            line = line.split()
            if "End" in line[0] and line[1] == block:
                print("Finished {} {}".format(block, name))
                return
            elif "Begin" in line[0] and line[1] not in self.SKIP_BLOCKS:
                # recursive call
                self._parse_block(fid, line[1], name=" ".join(line[2:]))
            elif "Default_Dir" in line[0]:
                directory = line[1]
            elif ".nc" in line[0]:
                file_path = self.file_path.parents[0] / directory / line[0]
                if directory:
                    data = self.data.setdefault(directory, {})
                else:
                    data = self.data
                data = data.setdefault(file_path.stem, {})
                data.update(parsers.parse_file("vlbi_netcdf", file_path=file_path).as_dict())
            else:
                data = self.data.setdefault(block, {})
                if name:
                    data = data.setdefault(name, {})
                data[line[0]] = " ".join(line[1:])
