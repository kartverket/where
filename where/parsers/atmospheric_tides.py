"""A parser for reading atmospheric tide coefficients

Description:
------------

Reads atmospheric tide coefficients. Gridded coefficients for the S1 and S2 atmospheric tides.

"""
# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_line import LineParser


@plugins.register
class AtmosphericTidesParser(LineParser):
    """A parser for reading gridded coefficents for atmospheric tides
    """

    def setup_parser(self):
        """Set up information needed for the parser

        This should return a dictionary with all parameters needed by np.genfromtxt to do the actual parsing.

        Returns:
            Dict:  Parameters needed by np.genfromtxt to parse the input file.
        """
        # No header or comments in the file
        return dict(
            names=[
                "lon",
                "lat",
                "A_d1_u",
                "B_d1_u",
                "A_d2_u",
                "B_d2_u",
                "A_d1_n",
                "B_d1_n",
                "A_d2_n",
                "B_d2_n",
                "A_d1_e",
                "B_d1_e",
                "A_d2_e",
                "B_d2_e",
            ]
        )
