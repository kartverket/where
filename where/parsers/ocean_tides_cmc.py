"""A parser for reading ocean tide center of mass coefficients

Description:
------------

Reads ocean tide center of mass coefficients.


Reference:
----------

http://holt.oso.chalmers.se/loading/cmc.html

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_line import LineParser


@plugins.register
class OceanTidesCmcParser(LineParser):
    """A parser for reading coefficents for ocean tides center of mass corrections
    """

    def setup_parser(self):
        """Set up information needed for the parser

        This should return a dictionary with all parameters needed by np.genfromtxt to do the actual parsing.

        Returns:
            Dict:  Parameters needed by np.genfromtxt to parse the input file.
        """
        return dict(
            skip_header=1,
            names=("tide", "model", "z_in", "z_cr", "x_in", "x_cr", "y_in", "y_cr"),
            dtype=("U3", "U42", "f8", "f8", "f8", "f8", "f8", "f8"),
        )

    def structure_data(self):
        """Make xyz-vectors of the coefficients
        """
        self.data["tide"] = self._array["tide"]
        self.data["model"] = self._array["model"]
        self.data["in_phase"] = np.stack((self._array["x_in"], self._array["y_in"], self._array["z_in"]), axis=1)
        self.data["cross_phase"] = np.stack((self._array["x_cr"], self._array["y_cr"], self._array["z_cr"]), axis=1)
