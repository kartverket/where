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
from where.parsers import parser
from where.lib import plugins
from where.ext import sofa

from where.parsers._parser import Parser

# Optional imports
from where.lib import optional

netCDF4 = optional.optional_import("netCDF4")


@plugins.register
class NetCDFParser(Parser):
    """A parser for reading ocean tidal loading coefficients from BLQ-files
    """

    SKIP_FIELDS_DEFAULT = [
        "Stub", "CreateTime", "CreatedBy", "Program", "Subroutine", "DataOrigin", "Session", "vgosDB_Version"
    ]

    def read_data(self):
        if True:
            self.SKIP_FIELDS = self.SKIP_FIELDS_DEFAULT
        else:
            self.SKIP_FIELDS = []

        data = netCDF4.Dataset(self.file_path)
        for key, variable in data.variables.items():
            if key in self.SKIP_FIELDS:
                continue

            self.data[key] = self._get_data(variable)

    def _get_data(self, variable):

        if variable.dtype == "S1":
            return netCDF4.chartostring(variable[:])
        else:
            return variable[:]
