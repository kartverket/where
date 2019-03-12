"""A parser for reading ocean tidal coefficients

Description:
------------

Parser for reading ocean tidal corrections to the gravity field C and S coefficients.

References:
-----------

    ftp://tai.bipm.org/iers/conv2010/chapter6/tidemodels/fes2004_Cnm-Snm.dat

"""

# Standard library imports
import re

# Midgard imports
from midgard.dev import plugins

from where.parsers._parser import Parser


@plugins.register
class OceanTidesCoeffParser(Parser):
    """A parser for reading ocean tidal coefficients
    """

    def read_data(self):
        with open(self.file_path) as fid:
            prog = re.compile("\d")
            for line in fid:
                # Search for line containing data
                if prog.match(line[2:3]):
                    self._parse_line(line, int(line[12:15]), int(line[16:19]))

    def _parse_line(self, line, n, m):
        self.data.setdefault("C+", {}).setdefault((n, m), {}).setdefault(float(line[0:7]), float(line[19:29]) * 1e-12)
        self.data.setdefault("C-", {}).setdefault((n, m), {}).setdefault(float(line[0:7]), float(line[29:39]) * 1e-12)
        self.data.setdefault("S+", {}).setdefault((n, m), {}).setdefault(float(line[0:7]), float(line[41:51]) * 1e-12)
        self.data.setdefault("S-", {}).setdefault((n, m), {}).setdefault(float(line[0:7]), float(line[51:61]) * 1e-12)
