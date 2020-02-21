"""A parser for reading VLBI data from vgosDb files

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
import netCDF4

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser import Parser


@plugins.register
class NetCDFParser(Parser):
    """A parser for reading netCDF files with VLBI data
    """

    SKIP_FIELDS_DEFAULT = [
        "Stub",
        "CreateTime",
        "CreatedBy",
        "Program",
        "Subroutine",
        "DataOrigin",
        "Session",
        "vgosDB_Version",
    ]

    def read_data(self):
        self.SKIP_FIELDS = self.SKIP_FIELDS_DEFAULT

        data = netCDF4.Dataset(self.file_path)
        for key, variable in data.variables.items():
            if key in self.SKIP_FIELDS:
                continue

            self.data[key] = self._get_data(variable)

    def _get_data(self, variable):
        variable.set_auto_mask(False)
        if variable.dtype == "S1":
            try:
                values = np.core.defchararray.strip(netCDF4.chartostring(variable[:]))
            except (UnicodeDecodeError, ValueError):
                # TODO: only happened with fields we are ignoring anyway so far
                values = np.core.defchararray.strip(netCDF4.chartostring(variable[:], encoding="bytes"))
        else:
            values = variable[:]

        if hasattr(variable, "REPEAT"):
            if values.ndim < 2:
                values = np.tile(values, variable.REPEAT)
            else:
                values = np.tile(values, variable.REPEAT).reshape(variable.REPEAT, -1)

        if "YMDHM" in variable.name:
            # Make sure year has 4 digits
            idx_add1900 = np.logical_and((values[:, 0] >= 50), (values[:, 0] < 100))
            idx_add2000 = values[:, 0] < 50
            if idx_add1900.any() or idx_add2000.any():
                values[:, 0][idx_add1900] += 1900
                values[:, 0][idx_add2000] += 2000

        return values
