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

    def _parse_line(self, line, n, m):
        limit = 0.05    # coefficient amplitude limit
        cp = float(line[19:29])
        sp = float(line[29:39])
        cm = float(line[41:51])
        sm = float(line[51:61])
        
        if ((abs(cp) < limit) and (abs(sp) < limit) and 
            (abs(cm) < limit) and (abs(sm) < limit)):
            return      # ignore entries with too low amplitude
        self.data.setdefault(float(line[0:7]), {}).setdefault("C+", {}).setdefault((n, m), cp * 1e-11)
        self.data.setdefault(float(line[0:7]), {}).setdefault("S+", {}).setdefault((n, m), sp * 1e-11)
        self.data.setdefault(float(line[0:7]), {}).setdefault("C-", {}).setdefault((n, m), cm * 1e-11)
        self.data.setdefault(float(line[0:7]), {}).setdefault("S-", {}).setdefault((n, m), sm * 1e-11)

