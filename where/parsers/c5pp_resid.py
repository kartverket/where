"""A parser for reading c5++ residuals

Description:
------------

Reads the orbit residuals that are output of the c5++ program

"""
# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_line import LineParser


@plugins.register
class C5ppResidParser(LineParser):
    """A parser for reading output files of c5pp
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
                "tech",
                "sat_name",
                "station",
                "_date",
                "_time",
                "time_scale",
                "calc",
                "residual",
                "good_bad",
                "_1",
                "_2",
                "_3",
                "_4",
                "_5",
                "_6",
                "_7",
                "_8",
                "_9",
                "_10",
                "_11",
                "_12",
                "_13",
            ],
            dtype="U3, U4, U4, U10, U16, U3, f8, f8, U1, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8",
        )
