"""A parser for reading atmospheric tide center of mass coefficients

Description:
------------

Reads atmospheric tide center of mass coefficients. 


Reference:
----------

http://geophy.uni.lu/ggfc-atmosphere/tide-loading-calculator.html




$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""
# Where imports
from where.parsers._parser_line import LineParser
from where.lib import plugins


@plugins.register
class AtmosphericTidesComParser(LineParser):
    """A parser for reading coefficents for atmospheric tides center of mass corrections
    """

    def setup_parser(self):
        """Set up information needed for the parser

        This should return a dictionary with all parameters needed by np.genfromtxt to do the actual parsing.

        Returns:
            Dict:  Parameters needed by np.genfromtxt to parse the input file.
        """
        return dict(names=True, dtype=("U2", "f8", "f8", "f8", "f8"))
