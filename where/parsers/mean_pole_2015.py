"""A parser for reading observed pole coordinates

Description:
------------


"""
# Midgard imports
from midgard.dev import plugins

# Where imports
from where.parsers._parser_line import LineParser


@plugins.register
class MeanPole2015(LineParser):
    """A parser for reading observed pole coordinates
    """

    def setup_parser(self):
        """Set up information needed for the parser

        This should return a dictionary with all parameters needed by np.genfromtxt to do the actual parsing.

        Returns:
            Dict:  Parameters needed by np.genfromtxt to parse the input file.
        """
        # No header or comments in the file
        return dict(names=["year", "x", "y"])
