"""A parser for reading ocean tidal coefficients

Description:
------------

Parser for reading ocean tidal corrections to the gravity field C and S coefficients.
from the S1 wave

References:
-----------

    ftp://tai.bipm.org/iers/conv2010/chapter6/tidemodels/S1.dat

"""

# Standard library imports
import re

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser import Parser


@plugins.register
class OceanTidesCoeffParser(Parser):
    """A parser for reading ocean tidal coefficients
    """

    def read_data(self):
        with open(self.file_path) as fid:
            prog = re.compile(r"\d")
            for line in fid:
                # Search for line containing data
                if prog.match(line[2:3]):
                    self._parse_line(line, int(line[12:15]), int(line[16:19]))

    # Doodson Darw  l   m    Csin+     Ccos+       Csin-     Ccos-       C+   eps+      C-   eps-
    # 164.556 S1    1   0 -0.023741  0.027064    0.000000  0.000000   0.0360 318.742 0.0000   0.000

    def _parse_line(self, line, n, m):
        self.data.setdefault(float(line[0:7]), {}).setdefault("C+", {}).setdefault((n, m), float(line[63:70]) * 1e-12)
        self.data.setdefault(float(line[0:7]), {}).setdefault("C-", {}).setdefault((n, m), float(line[79:85]) * 1e-12)
