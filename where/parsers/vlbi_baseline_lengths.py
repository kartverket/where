"""A parser for reading baseline lengths

Description:
------------

The baseline lengths file is produced by the postprosessor vlbi_baseline_length

"""

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_line import LineParser


@plugins.register
class BaselineLengths(LineParser):
    """A parser for reading baseline lengths
    """

    def setup_parser(self):
        """Set up information needed for the parser

        This should return a dictionary with all parameters needed by np.genfromtxt to do the actual parsing.

        Returns:
            Dict:  Parameters needed by np.genfromtxt to parse the input file.
        """
        return dict(names=["baseline", "num_obs", "length", "ferr"],
                    dtype=["U20", "i8", "f8", "f8"])
    
    def structure_data(self):
        for item in self._array:
            self.data[item["baseline"]] = dict(length=item["length"], ferr=item["ferr"], num_obs=item["num_obs"])
