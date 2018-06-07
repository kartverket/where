"""A parser for reading data from ITRF files in SNX format

Description:
------------

Reads station positions and velocities from ITRF files in SNX format.




"""

# Standard library imports
from datetime import datetime, timedelta
import itertools

# Where imports
from where.lib import plugins
from where.lib import config
from where.parsers._parser_sinex import SinexParser, SinexBlock, SinexField, parsing_factory


@plugins.register
class VlbiSnxEstimateParser(SinexParser):
    """A parser for reading data from ITRF files in SNX format
    """

    def setup_parser(self):
        return (self.solution_estimate, self.solution_estimates, self.solution_apriori)

    @property
    def solution_estimates(self):
        """Estimated parameters.

        Some analysis centers call this block SOLUTION/ESTIMATES instead of SOLUTION/ESTIMATE. The content of the block
        is otherwise according the the format specification.

        Example:
            *INDEX TYPE__ CODE PT SOLN _REF_EPOCH__ UNIT S __ESTIMATED VALUE____ _STD_DEV___
                 1 STAX   7207  A    1 10:001:00000 m    2 -.240960109141758E+07 0.12784E-02
                      1111111111222222222233333333334444444444555555555566666666667777777777
            01234567890123456789012345678901234567890123456789012345678901234567890123456789
        """
        return SinexBlock(
            marker="SOLUTION/ESTIMATES",
            fields=(
                SinexField("param_idx", 1, "i8"),
                SinexField("param_name", 7, "U6"),
                SinexField("site_code", 14, "U4"),
                SinexField("point_code", 19, "U2"),
                SinexField("soln", 22, "U4"),
                SinexField("ref_epoch", 27, "O", "epoch"),
                SinexField("unit", 40, "U4"),
                SinexField("constraint", 45, "U1"),
                SinexField("estimate", 47, "f8", "exponent"),
                SinexField("estimate_std", 69, "f8", "exponent"),
            ),
            parser=self.parse_solution_estimates,
        )

    def parse_solution_estimates(self, data):
        self.data["SOLUTION/ESTIMATE"] = self.parse_solution_estimate(data)
